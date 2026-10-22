#!/usr/bin/env python3
"""cnmf_fit.py — cNMF arm for the MR-2 benchmark (cnmf 1.7.1 API).

Reads a donor x gene counts matrix (mtx or mtx.gz), runs cNMF over a K grid,
selects K by the consensus-stability criterion registered in
PREREG_mr2-cnmf-v1.md (largest K whose mean stability is within 0.05 of the
grid maximum; ties -> smaller K), and writes usage / spectra / chosen-K.

cnmf 1.7.1 exposes stability as the consensus *silhouette* (higher = stabler),
not the older local_density. The selection rule is applied to silhouette.

Usage: python3 cnmf_fit.py <counts.mtx[.gz]> <out_prefix> [--kgrid 5,10,20,30]
Outputs: <out_prefix>.usage.tsv  (donors x K)
         <out_prefix>.spectra.tsv (genes x K)
         <out_prefix>.chosen_k.txt
         <out_prefix>.stability.csv
"""
import sys
import os
import gzip
import tempfile

import numpy as np
import pandas as pd
import scipy.io as scio
import scipy.sparse as sp


def read_counts(path):
    op = gzip.open if path.endswith(".gz") else open
    M = scio.mmread(op(path))
    return M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)


def main(counts_path, out_prefix, kgrid=(5, 10, 20, 30)):
    from cnmf import cNMF
    X = read_counts(counts_path)  # donors x genes, counts
    donors = [f"d{i}" for i in range(X.shape[0])]
    genes = [f"g{i}" for i in range(X.shape[1])]
    with tempfile.TemporaryDirectory() as td:
        import anndata as ad
        obs = pd.DataFrame(index=donors)
        var = pd.DataFrame(index=genes)
        a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
        counts_fn = os.path.join(td, "counts.h5ad")
        a.write_h5ad(counts_fn)

        obj = cNMF(output_dir=td, name="cnmf")
        obj.prepare(counts_fn=counts_fn, components=list(kgrid), n_iter=100,
                    seed=14, num_highvar_genes=min(2000, X.shape[1]))
        obj.factorize(worker_i=0, total_workers=1)
        obj.combine()

        # per-K stability (silhouette) + prediction error, no figures
        # (consensus reads the prepare()-saved normalized counts when norm_counts=None)
        stats = []
        for k in sorted(set(int(x) for x in kgrid)):
            st = obj.consensus(k, skip_density_and_return_after_stats=True,
                               show_clustering=False, close_clustergram_fig=True).stats
            stats.append({"k": k, "silhouette": float(st["silhouette"]),
                          "prediction_error": float(st["prediction_error"])})
        stats = pd.DataFrame(stats).sort_values("k").reset_index(drop=True)

        # registered selection: largest K with stability >= max - 0.05; ties -> smaller K
        smax = stats["silhouette"].max()
        eligible = stats[stats["silhouette"] >= smax - 0.05]
        chosen = int(eligible["k"].max())

        obj.consensus(k=chosen, density_threshold=0.2, show_clustering=False,
                      close_clustergram_fig=True)
        usage, spectra_scores, spectra_tpm, top_genes = obj.load_results(
            chosen, 0.2, norm_usage=False)

        # load_results returns usage (donors x K) and spectra_tpm (genes x K,
        # already transposed). Re-index to our donor/gene labels and standardize
        # column names to P1..PK.
        usage = pd.DataFrame(np.asarray(usage), index=donors,
                             columns=[f"P{i+1}" for i in range(usage.shape[1])])
        spectra = pd.DataFrame(np.asarray(spectra_tpm), index=genes,
                               columns=[f"P{i+1}" for i in range(spectra_tpm.shape[1])])

        usage.to_csv(out_prefix + ".usage.tsv", sep="\t")
        spectra.to_csv(out_prefix + ".spectra.tsv", sep="\t")
        with open(out_prefix + ".chosen_k.txt", "w") as fh:
            fh.write(f"{chosen}\n")
        stats.to_csv(out_prefix + ".stability.csv", index=False)
        print(f"chosen K={chosen}; silhouette: " +
              ", ".join(f"k{int(r.k)}={r.silhouette:.3f}" for r in stats.itertuples()))


if __name__ == "__main__":
    a = sys.argv[1:]
    kg = (5, 10, 20, 30)
    if "--kgrid" in a:
        i = a.index("--kgrid")
        kg = tuple(int(x) for x in a[i + 1].split(","))
        a = a[:i] + a[i + 2:]
    main(a[0], a[1], kg)
