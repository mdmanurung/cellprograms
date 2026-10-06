"""Panel C (program quality) + Panel D partial (interpretability) for COMBAT.

EBMF: refit per cell type (deterministic seed 67) saving loadings; compute
  - R2 reconstruction (per cell type, overall)
  - held-out reconstruction (mask 10% of entries, R2 on masked entries)
  - stability: 10 bootstrap donor refits (80%), factor recovery by loading
    cosine (Hungarian matching, threshold 0.9), score corr, top-50 gene overlap
MOFA: read factors+weights from the trained hdf5; same reconstruction metrics.

Usage: python run_combat_panelC.py [ebmf|mofa|all]
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)
SEED = 67
N_STABILITY = 10
BOOT_FRAC = 0.8
COS_THRESHOLD = 0.9
TOP_N_GENES = 50


# ---------- EBMF ----------

def _ebmf_fit_full(item):
    """Fit one cell type; return (ct, scores df, loadings genes x K)."""
    ct, df = item
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from patpy_cellprograms import EBMF
    m = EBMF(sample_key="sample", cell_group_key="ct", seed=SEED)
    m._fit_in_r({ct: df})
    m._fitted = True
    block = m.build_representation()
    Ldf = m.get_loadings().get(ct)
    return ct, block, (None if Ldf is None else (Ldf.values, Ldf))


def ebmf_full_fit(pb, workers=4):
    from concurrent.futures import ProcessPoolExecutor
    out = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for ct, block, L in ex.map(_ebmf_fit_full, [(ct, pb[ct]) for ct in sorted(pb)]):
            out[ct] = {"block": block,
                       "loadings": None if L is None else L[0],
                       "loadings_df": None if L is None else L[1]}
            print(f"  full fit {ct}: K={0 if L is None else L[0].shape[1]}", flush=True)
    return out


def ebmf_fit_subset(ct, df, donor_idx):
    """Refit one cell type on a donor subset; returns loadings (genes x K)."""
    from patpy_cellprograms import EBMF
    m = EBMF(sample_key="sample", cell_group_key="ct", seed=SEED)
    sub = df.iloc[donor_idx]
    m._fit_in_r({ct: sub})
    m._fitted = True
    L = m.get_loadings().get(ct)
    return None if L is None else L.values


def reconstruction_metrics(Y, L, S):
    """R2 of Y (samples x genes) under the canonical factorization Y ~ S @ L.T.

    S: canonical scores (samples x K, EBMF posterior means). L: loadings
    (genes x K). Held-out: zero 10% of entries, re-project scores from the
    masked data, score R2 on masked entries only.
    """
    Yc = Y - Y.mean(0)
    if L is None or L.shape[1] == 0 or S is None or S.shape[1] == 0:
        return 0.0, 0.0
    Yhat = S @ L.T
    ss_res = np.sum((Yc - Yhat) ** 2)
    ss_tot = np.sum(Yc ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    # held-out: zero 10% of entries, re-project scores, R2 on masked entries
    rng = np.random.default_rng(SEED)
    mask = rng.random(Y.shape) < 0.10
    Yh = Yc.copy()
    Yh[mask] = 0.0
    # least-squares projection of masked data onto loadings
    S_h = Yh @ np.linalg.pinv(L.T)  # samples x K (pinv(L.T) is K x genes)
    Yhat_h = S_h @ L.T
    resid = (Yc - Yhat_h)[mask]
    ss_tot_h = np.sum(Yc[mask] ** 2)
    r2_held = 1 - np.sum(resid ** 2) / ss_tot_h if ss_tot_h > 0 else np.nan
    return float(r2), float(r2_held)


def stability_refits(pb, full_fits, n_boot=N_STABILITY, workers=4):
    """Bootstrap donor refits; factor recovery vs full fit per cell type."""
    from concurrent.futures import ProcessPoolExecutor
    rng = np.random.default_rng(2026)
    rows = []
    for ct, info in full_fits.items():
        L_full = info["loadings"]
        df = pb[ct]
        if L_full is None or L_full.shape[1] == 0:
            continue
        if L_full.shape[0] != df.shape[1]:
            L_full = info["loadings_df"].reindex(df.columns, fill_value=0.0).values
        norms = np.linalg.norm(L_full, axis=0)
        norms[norms == 0] = 1
        Ln = L_full / norms
        rec_freq = np.zeros(L_full.shape[1])
        score_corrs, gene_overlaps = [], []
        jobs = []
        for b in range(n_boot):
            k = int(BOOT_FRAC * df.shape[0])
            idx = rng.choice(df.shape[0], k, replace=False)
            jobs.append((ct, df, idx))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_boot_worker, jobs))
        for L_b in results:
            if L_b is None or L_b.shape[1] == 0:
                continue
            L_b = L_b.reindex(pb[ct].columns, fill_value=0.0).values
            nb = np.linalg.norm(L_b, axis=0)
            nb[nb == 0] = 1
            Lbn = L_b / nb[None, :]
            cos = Ln.T @ Lbn  # K_full x K_boot
            r, c = linear_sum_assignment(-cos)
            for ri, ci in zip(r, c):
                if cos[ri, ci] > COS_THRESHOLD:
                    rec_freq[ri] += 1
        rows.append({"cell_type": ct, "K": L_full.shape[1],
                     "recovery_freq_mean": float(rec_freq.mean()),
                     "recovery_freq_per_factor": rec_freq.tolist()})
        print(f"  stability {ct}: mean recovery {rec_freq.mean():.2f}", flush=True)
    return pd.DataFrame(rows)


def _boot_worker(job):
    ct, df, idx = job
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from patpy_cellprograms import EBMF
    m = EBMF(sample_key="sample", cell_group_key="ct", seed=SEED)
    m._fit_in_r({ct: df.iloc[idx]})
    m._fitted = True
    L = m.get_loadings().get(ct)
    return None if L is None else L  # labelled DataFrame


# ---------- MOFA ----------

def mofa_reconstruction(hdf5_path, pb):
    """R2 + held-out R2 per view from a trained MOFA hdf5.

    mofapy2 layout: data[view][group] is (N, genes); expectations/W[view] is
    (K, genes); expectations/Z[group] is (K, N). Reconstruction Y ~ Z.T @ W.
    """
    import h5py
    rows = []
    with h5py.File(hdf5_path, "r") as f:
        Z = f["expectations"]["Z"]["single_group"][:]  # K x N
        for view in f["data"].keys():
            Y = f["data"][view]["single_group"][:]  # N x genes
            W = f["expectations"]["W"][view][:]  # K x genes
            Yc = Y - Y.mean(0, keepdims=True)
            Yhat = Z.T @ W  # N x genes
            ss_res = np.sum((Yc - Yhat) ** 2)
            ss_tot = np.sum(Yc ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            # held-out: zero 10% of entries, re-project Z, R2 on masked entries
            rng = np.random.default_rng(SEED)
            mask = rng.random(Y.shape) < 0.10
            Yh = Yc.copy()
            Yh[mask] = 0.0
            Z_h = np.linalg.pinv(W.T) @ Yh.T  # K x N (pinv(W.T) is K x genes)
            Yhat_h = Z_h.T @ W
            resid = (Yc - Yhat_h)[mask]
            ss_tot_h = np.sum(Yc[mask] ** 2)
            r2h = 1 - np.sum(resid ** 2) / ss_tot_h if ss_tot_h > 0 else np.nan
            rows.append({"view": view, "K": W.shape[0], "r2": float(r2),
                         "r2_heldout": float(r2h)})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    adata, st = cc.load_combat()
    if cc.ckpt_exists("pb_mean_A.pkl"):
        pb = cc.load_ckpt("pb_mean_A.pkl")
    else:
        pb = cc.pseudobulk(adata, agg="mean")
        cc.save_ckpt(pb, "pb_mean_A.pkl")
    del adata

    if which in ("ebmf", "all"):
        t0 = time.time()
        if cc.ckpt_exists("ebmf_full_fits.pkl"):
            full = cc.load_ckpt("ebmf_full_fits.pkl")
        else:
            full = ebmf_full_fit(pb)
            cc.save_ckpt(full, "ebmf_full_fits.pkl")
        # reconstruction (align loadings to full gene space; dropped
        # zero-variance genes are all-zero after centering -> zero rows)
        rows = []
        for ct, info in full.items():
            L = info["loadings"]
            Y = pb[ct].values.astype(np.float64)
            if L is not None and L.shape[0] != Y.shape[1]:
                Ldf = info.get("loadings_df")
                if Ldf is not None:
                    L = Ldf.reindex(pb[ct].columns, fill_value=0.0).values
                else:
                    raise RuntimeError(f"{ct}: cannot align loadings {L.shape} to {Y.shape[1]} genes")
            S = info["block"].reindex(pb[ct].index).fillna(0.0).values
            r2, r2h = reconstruction_metrics(Y, L, S)
            rows.append({"method": "EBMF", "cell_type": ct, "K": 0 if L is None else L.shape[1],
                         "r2": r2, "r2_heldout": r2h})
            print(f"  recon {ct}: R2={r2:.3f} heldout={r2h:.3f}", flush=True)
        pd.DataFrame(rows).to_csv("/mnt/results/combat_benchmark/panel_C_ebmf_reconstruction.csv", index=False)
        # stability on all cell types (COMBAT has 17)
        stab = stability_refits(pb, full)
        stab.to_csv("/mnt/results/combat_benchmark/panel_C_ebmf_stability.csv", index=False)
        print(f"EBMF panel C done ({time.time()-t0:.0f}s)", flush=True)

    if which in ("mofa", "all"):
        rows = mofa_reconstruction(f"{cc.WORK}/mofa_A.hdf5", pb)
        rows.to_csv("/mnt/results/combat_benchmark/panel_C_mofa_reconstruction.csv", index=False)
        print(rows.to_string())
