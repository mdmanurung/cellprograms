"""M3 S5 anchor: cNMF on the same synthetic matrices (pre-registered 04 §1).

cNMF (python package) with density-silhouette K selection over 2..8, 50 restarts
(anchor role, compute-motivated; E's 100-restart amendment targeted K-selection
fidelity in the detection pilot - the deviation is recorded). Metrics vs planted
truth: Hungarian |cor| of usage vs Z_true on train cells; top-50 Jaccard of spectra
vs marker blocks (max possible = 50/len(markers)).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from cnmf import cNMF

from m3_metrics import hungarian_match, top_jaccard
from m3_synth import generate

warnings.filterwarnings("ignore")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--restarts", type=int, default=50)
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--tmp", default="tmp_cnmf")
    args = ap.parse_args()

    g = generate(args.rep, args.rate, args.beta)
    train_idx = np.where(~g["test_mask"])[0]
    counts = pd.DataFrame(g["counts"][train_idx].toarray(),
                          index=[f"c{i}" for i in range(len(train_idx))],
                          columns=[f"G{j}" for j in range(g["counts"].shape[1])])
    marker_sets = [set(m) for m in g["markers"]]

    os.makedirs(args.tmp, exist_ok=True)
    name = f"m3_rate{args.rate}_rep{args.rep}"
    h5 = os.path.join(args.tmp, name + ".h5ad")
    ad.AnnData(X=counts.values.astype(np.float32),
               obs=pd.DataFrame(index=counts.index),
               var=pd.DataFrame(index=counts.columns)).write_h5ad(h5)
    obj = cNMF(output_dir=args.tmp, name=name)
    t0 = time.time()
    obj.prepare(counts_fn=h5, components=list(range(2, 9)), seed=args.rep,
                n_iter=args.restarts, num_highvar_genes=1500)
    obj.factorize()
    obj.combine(components=list(range(2, 9)))
    t_fit = time.time() - t0

    t0 = time.time()
    sil = {}
    for k in range(2, 9):
        try:
            st = obj.consensus(k=k, skip_density_and_return_after_stats=True,
                               show_clustering=False, close_clustergram_fig=True).stats
            v = st.iloc[0]["silhouette"] if isinstance(st, pd.DataFrame) else st["silhouette"]
            sil[k] = float(v)
        except Exception:
            sil[k] = float("nan")
    k_hat = max(range(2, 9), key=lambda k: (sil[k] if sil[k] == sil[k] else -1, -k))
    obj.consensus(k=k_hat, show_clustering=False, close_clustergram_fig=True)
    usage, spectra, _, _ = obj.load_results(K=k_hat, density_threshold=0.5)
    t_sel = time.time() - t0

    Z_train = g["Z_true"][train_idx]
    usage_al = usage.loc[counts.index].values
    matches = hungarian_match(usage_al, Z_train)
    for m in matches:
        m["top50_jaccard"] = top_jaccard(spectra.values, m["est_col"],
                                         marker_sets[m["truth_col"]])
    row = dict(arm="S5", rate=args.rate, rep=args.rep, beta=args.beta,
               n_train=len(train_idx), k_hat=int(k_hat),
               silhouette={str(k): round(v, 4) for k, v in sil.items()},
               t_fit_s=round(t_fit, 2), t_sel_s=round(t_sel, 2), matches=matches,
               zero_rate=float((counts.values == 0).mean()), restarts=args.restarts)
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, f"S5_rate{args.rate}_rep{args.rep:02d}.json"), "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: v for k, v in row.items() if k not in ("matches", "silhouette")}))
    print("matches:", [(m["truth_col"], round(m["match_cor"], 3),
                         round(m["top50_jaccard"], 2)) for m in matches])
    shutil.rmtree(os.path.join(args.tmp, name), ignore_errors=True)
    os.remove(h5)


if __name__ == "__main__":
    main()
