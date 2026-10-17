# Workstream E: cNMF arm of the ebmf-vs-cNMF audit pilot (v3).
# v3 changes (protocol amendment, documented in TASKS.md before unblinded
# analysis): restarts raised to the cNMF default 100 for BOTH arms (the
# 10/15-restart deviation was compute-motivated; measured cost 0.05 s/task
# made it unnecessary, and at 10 restarts cNMF's built-in silhouette
# misselected K on smoke rep 1). K selection: cNMF's own density-silhouette
# over K=2..8 (undefined at K=1; cNMF has no K=0 concept — noise-scenario
# scope limitation reported separately). Supplementary restart-stability
# (mean pairwise greedy max-matched |cor| across restart spectra) recorded
# but not used for selection.
import argparse
import itertools
import os
import shutil
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz

warnings.filterwarnings("ignore")

BLOCK = 30


def pearson(a, b):
    return float(np.corrcoef(np.asarray(a, dtype=float), np.asarray(b, dtype=float))[0, 1])


def matched_mean_cor(A, B):
    C = np.abs(np.corrcoef(A, B)[:A.shape[0], A.shape[0]:])
    C = np.nan_to_num(C, nan=0.0)
    total, n = 0.0, A.shape[0]
    for _ in range(n):
        i, j = np.unravel_index(np.argmax(C), C.shape)
        total += C[i, j]
        C = np.delete(C, i, axis=0)
        C = np.delete(C, j, axis=1)
    return total / n


def restart_stability(obj, ks, n_iter):
    out = {}
    for K in ks:
        mats = [load_df_from_npz(obj.paths["iter_spectra"] % (K, r)).values
                for r in range(n_iter)]
        vals = [matched_mean_cor(mats[a], mats[b])
                for a, b in itertools.combinations(range(n_iter), 2)]
        out[K] = float(np.mean(vals))
    return out


def _stat(st, key):
    try:
        if isinstance(st, pd.DataFrame):
            v = st.iloc[0][key]
        elif isinstance(st, dict):
            v = st[key]
        else:
            v = getattr(st, key)
        return float(v)
    except Exception:
        return float("nan")


def cnmf_silhouette(obj, k):
    st = obj.consensus(k=k, skip_density_and_return_after_stats=True,
                       show_clustering=False, close_clustergram_fig=True).stats
    return _stat(st, "silhouette"), _stat(st, "prediction_error")


def fit_cnmf(counts, out_dir, name, components, seed, n_restarts, num_hvg):
    h5 = os.path.join(out_dir, name + ".h5ad")
    ad.AnnData(
        X=counts.values.astype(np.float32),
        obs=pd.DataFrame(index=counts.index),
        var=pd.DataFrame(index=counts.columns),
    ).write_h5ad(h5)
    obj = cNMF(output_dir=out_dir, name=name)
    t0 = time.time()
    obj.prepare(counts_fn=h5, components=components, seed=seed, n_iter=n_restarts,
                num_highvar_genes=num_hvg)
    obj.factorize()
    obj.combine(components=components)
    t_fit = time.time() - t0
    return obj, t_fit


def cleanup(out_dir, name):
    d = os.path.join(out_dir, name)
    if os.path.isdir(d):
        shutil.rmtree(d, ignore_errors=True)
    f = os.path.join(out_dir, name + ".h5ad")
    if os.path.exists(f):
        os.remove(f)


def greedy_match(cors):
    matches = []
    used_c = set()
    order = sorted(((cors[t, c], t, c) for t in range(cors.shape[0])
                    for c in range(cors.shape[1])), reverse=True)
    taken_t = set()
    for v, t, c in order:
        if t in taken_t or c in used_c:
            continue
        matches.append((t, c, v))
        taken_t.add(t)
        used_c.add(c)
    return matches


def run_recovery(sc, rep, k_true, data_dir, out_dir):
    tag = "recovery_%s_%02d" % (sc, rep)
    counts = pd.read_csv(os.path.join(data_dir, tag + "_counts.csv"), index_col=0)
    ks = list(range(1, 9))
    obj, t_fit = fit_cnmf(counts, out_dir, tag, ks, seed=rep, n_restarts=100,
                          num_hvg=600)
    t0 = time.time()
    sil = {}
    for k in ks[1:]:
        try:
            sil[k] = cnmf_silhouette(obj, k)
        except Exception:
            sil[k] = (float("nan"), float("nan"))
    rstab = restart_stability(obj, ks, n_iter=100)
    t_sel = time.time() - t0
    k_hat = max(ks[1:], key=lambda k: (sil[k][0], -k))

    obj.consensus(k=k_hat, show_clustering=False, close_clustergram_fig=True)
    usage, spectra, _, _ = obj.load_results(K=k_hat, density_threshold=0.5)

    row = {"scenario": sc, "rep": rep, "k_true": k_true, "k_hat": k_hat, "t_fit_s": t_fit,
           "t_sel_s": t_sel, "silhouette_khat": sil[k_hat][0],
           "max_silhouette": max(v[0] for v in sil.values()),
           "restart_stab_khat": rstab[k_hat]}
    if k_true > 0:
        S = pd.read_csv(os.path.join(data_dir, tag + "_scores.csv"), index_col=0)
        cors = np.zeros((k_true, usage.shape[1]))
        for k in range(k_true):
            for j in range(usage.shape[1]):
                cors[k, j] = abs(pearson(usage.iloc[:, j], S.iloc[:, k]))
        matches = greedy_match(cors)
        genes = counts.columns.tolist()
        for t, c, v in matches:
            block = set(genes[BLOCK * t:BLOCK * (t + 1)])
            top30 = set(spectra.iloc[:, c].nlargest(30).index)
            row[f"score_cor_P{t + 1}"] = v
            row[f"jaccard_P{t + 1}"] = len(block & top30) / len(block | top30)
    cleanup(out_dir, tag)
    return row


def run_detection(rep, rho, data_dir, out_dir, rate=0.0):
    if rate > 0:
        tag = "det_drop%g_rho%g_%02d" % (rate, rho, rep)
    else:
        tag = "det_rho%g_%02d" % (rho, rep)
    usages = {}
    t_fit = 0.0
    for ct in ["A", "B"]:
        counts = pd.read_csv(os.path.join(data_dir, "%s_%s_counts.csv" % (tag, ct)),
                             index_col=0)
        name = tag + "_" + ct
        obj, tf = fit_cnmf(counts, out_dir, name, [10], seed=rep, n_restarts=100,
                           num_hvg=600)
        t_fit += tf
        obj.consensus(k=10, show_clustering=False, close_clustergram_fig=True)
        usage, _, _, _ = obj.load_results(K=10, density_threshold=0.5)
        usages[ct] = usage
        cleanup(out_dir, name)

    row = {"rep": rep, "rate": rate, "rho": rho, "t_fit_s": t_fit}
    S = {}
    for ct in ["A", "B"]:
        # dropout data shares the clean data's seed stream -> identical truth
        # scores live in the clean tags' score CSVs
        scores_tag = ("det_rho%g_%02d" % (rho, rep)) if rate > 0 else tag
        S[ct] = pd.read_csv(os.path.join(data_dir, "%s_scores_%s.csv" % (scores_tag, ct)),
                            index_col=0)
    cors_A = [abs(pearson(usages["A"].iloc[:, j], S["A"].iloc[:, 0])) for j in range(10)]
    cors_B = [abs(pearson(usages["B"].iloc[:, j], S["B"].iloc[:, 0])) for j in range(10)]
    jA, jB = int(np.argmax(cors_A)), int(np.argmax(cors_B))
    row["match_A"] = cors_A[jA]
    row["match_B"] = cors_B[jB]
    row["detected"] = int(cors_A[jA] > 0.5 and cors_B[jB] > 0.5)
    if row["detected"]:
        row["est_cor"] = pearson(usages["A"].iloc[:, jA], usages["B"].iloc[:, jB])
    C = np.corrcoef(np.hstack([usages["A"].values, usages["B"].values]))[:10, 10:]
    n_pairs = int((np.abs(C) > 0.5).sum())
    row["n_false_pairs"] = max(n_pairs - (1 if row["detected"] else 0), 0)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=["recovery", "detection"], required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--scenarios", default="k3,noise",
                    help="recovery scenarios: k3, noise, drop0.1, drop0.3")
    ap.add_argument("--drop-rate", type=float, default=0.0,
                    help="detection dropout rate (0 = clean detection data)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    if args.arm == "recovery":
        spec = {"k3": (range(1, 11), 3), "noise": (range(11, 21), 0),
                "drop0.1": (range(1, 11), 3), "drop0.3": (range(1, 11), 3)}
        for sc in args.scenarios.split(","):
            reps, k_true = spec[sc]
            for rep in reps:
                rows.append(run_recovery(sc, rep, k_true, args.data_dir,
                                         args.out_dir))
                pd.DataFrame(rows).to_csv(args.out_csv, index=False)
    else:
        combos = [(rep, rho) for rep in range(1, 21)
                  for rho in [0, 0.3, 0.45, 0.6]]
        combos = [c for i, c in enumerate(combos) if i % args.nshards == args.shard]
        for rep, rho in combos:
            rows.append(run_detection(rep, rho, args.data_dir, args.out_dir,
                                      rate=args.drop_rate))
            pd.DataFrame(rows).to_csv(args.out_csv, index=False)
    print("done: %d rows -> %s" % (len(rows), args.out_csv))


if __name__ == "__main__":
    main()
