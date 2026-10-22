#!/usr/bin/env python3
"""cnmf_fit.py -- cNMF factorization for one cell type (MR-2 arm).

Reads a donor x gene counts matrix (mtx.gz or plain mtx), runs cNMF over the
registered K grid {5,10,20,30}, selects K by the registered stability rule
operationalized as: stability_k = 1 - local_density_k / max(local_density)
(cNMF's k-selection metric, rescaled to [0,1), higher = more stable); pick
the LARGEST k with stability_k >= max(stability) - 0.05 (ties -> smaller k).

Writes: <out_prefix>.usage.tsv (donors x K), <out_prefix>.spectra.tsv
(genes x K), <out_prefix>.chosen_k.txt, <out_prefix>.stability.csv.

Usage: python3 cnmf_fit.py <counts.mtx[.gz]> <out_prefix> [--kgrid 5,10,20,30]
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import io as scio


def read_counts(path):
    op = gzip.open if path.endswith(".gz") else open
    M = scio.mmread(path) if not path.endswith(".gz") else scio.mmread(op(path))
    M = M.tocsr()
    return M


def main(counts_path, out_prefix, kgrid=(5, 10, 20, 30)):
    from cnmf import cNMF
    X = read_counts(counts_path)  # donors x genes, counts
    donors = [f"d{i}" for i in range(X.shape[0])]
    genes = [f"g{i}" for i in range(X.shape[1])]
    with tempfile.TemporaryDirectory() as td:
        import anndata as ad
        obs = pd.DataFrame(index=donors)
        var = pd.DataFrame(index=genes)
        ad.X = None
        a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
        counts_fn = os.path.join(td, "counts.h5ad")
        a.write_h5ad(counts_fn)
        obj = cNMF(td, n_workers=1)
        obj.prepare(counts_fn=counts_fn, t=1, seed=14,
                    k_search_space=",".join(map(str, kgrid)), num_iter=100)
        obj.factorize()
        obj.k_selection_plot(close_figures=False)
        ld = obj.local_density.sort_index()  # k -> local density (lower = stabler)
        stab = 1.0 - ld / ld.max()
        chosen = None
        for k in sorted(ld.index, reverse=True):  # largest first
            if stab.loc[k] >= stab.max() - 0.05:
                chosen = int(k)
                break
        obj.consensus(k=chosen, density_threshold=0.2, local_density_threshold=0.05)
        usage, spectra = obj.get_results(k=chosen)
        usage.index = donors
        spectra.index = genes
        usage.to_csv(out_prefix + ".usage.tsv", sep="\t")
        spectra.to_csv(out_prefix + ".spectra.tsv", sep="\t")
        with open(out_prefix + ".chosen_k.txt", "w") as fh:
            fh.write(f"{chosen}\n")
        pd.DataFrame({"k": ld.index, "local_density": ld.values,
                      "stability": stab.values}).to_csv(out_prefix + ".stability.csv",
                                                        index=False)
        print(f"chosen K={chosen}; stability: " +
              ", ".join(f"k{k}={stab.loc[k]:.3f}" for k in ld.index))


if __name__ == "__main__":
    a = sys.argv[1:]
    kg = (5, 10, 20, 30)
    if "--kgrid" in a:
        i = a.index("--kgrid")
        kg = tuple(int(x) for x in a[i + 1].split(","))
        a = a[:i] + a[i + 2:]
    main(a[0], a[1], kg)
