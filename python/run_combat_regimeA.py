"""Regime A fits: EBMF (parallel per-cell-type via adapter), MOFAcellulaR matched-input.

Usage: python run_combat_regimeA.py [ebmf|mofa|all]
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)


def _ebmf_worker(item):
    """Subprocess: fit one cell type through the adapter's proven R path."""
    ct, df = item
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from patpy_cellprograms import EBMF

    m = EBMF(sample_key="sample", cell_group_key="ct", seed=67)
    m._fit_in_r({ct: df})  # single-cell-type fit, exact adapter code path
    m._fitted = True
    block = m.build_representation()  # z-scaled scores block for this ct
    return ct, block


def fit_ebmf(pb, workers=4):
    from concurrent.futures import ProcessPoolExecutor

    cts = sorted(pb.keys())
    reps = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for ct, block in ex.map(_ebmf_worker, [(ct, pb[ct]) for ct in cts]):
            reps.append(block)
            print(f"  EBMF done: {ct} ({block.shape[1]} factors, {time.time()-t0:.0f}s)", flush=True)
    rep = pd.concat(reps, axis=1).fillna(0.0)
    rep.columns = rep.columns.astype(str)
    return rep


def fit_mofa(pb, n_factors=10, seed=67, tag="A", hvg=False):
    """MOFAcellulaR via standalone Rscript (mofa_fit.R).

    tag="A": matched input — log-normalized mean pseudobulk as logcounts assay.
    tag="B"/"B2": raw-count sum pseudobulk -> filt_profiles -> tmm_trns
    (-> scran top-2000 HVG per view if hvg) -> MOFA2.
    """
    import subprocess, tempfile

    frames, donor, ctype = [], [], []
    for ct, df in pb.items():
        d = df.T.copy()  # genes x samples
        d.columns = [f"{ct}__{s}" for s in df.index]
        frames.append(d)
        donor += list(df.index)
        ctype += [ct] * df.shape[0]
    mat = pd.concat(frames, axis=1)
    coldata = pd.DataFrame({"donor_id": donor, "cell_type": ctype, "cell_counts": 100},
                           index=mat.columns)
    tmp = tempfile.mkdtemp(dir=cc.WORK)
    mat_path, cd_path = f"{tmp}/mat.rds", f"{tmp}/coldata.rds"
    out_path = f"{tmp}/factors_{tag}.csv"
    save_rds(mat, mat_path)
    save_rds(coldata, cd_path)
    cmd = ["Rscript", "mofa_fit.R", mat_path, cd_path, tag, str(n_factors),
           str(seed), "1" if hvg else "0", out_path]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True,
                         env={**os.environ, "R_HOME": os.popen("R RHOME").read().strip()})
    print(res.stdout[-2000:], flush=True)
    if res.returncode != 0:
        raise RuntimeError(f"MOFA fit {tag} failed:\n{res.stderr[-3000:]}")
    print(f"MOFA-{tag}: {time.time()-t0:.0f}s", flush=True)
    fac = pd.read_csv(out_path, index_col=0)
    return fac


def save_rds(df, path):
    """Write a DataFrame to RDS via a minimal Rscript call."""
    import subprocess, tempfile
    tmp_csv = path.replace(".rds", ".csv")
    df.to_csv(tmp_csv)
    script = (f'x <- read.csv("{tmp_csv}", row.names=1, check.names=FALSE); '
              f'saveRDS(x, "{path}")')
    subprocess.run(["Rscript", "-e", script], check=True, capture_output=True,
                   env={**os.environ, "R_HOME": os.popen("R RHOME").read().strip()})


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    adata, st = cc.load_combat()
    print(f"loaded: {adata.shape}, {len(st)} samples ({time.time()-t0:.0f}s)", flush=True)
    if cc.ckpt_exists("pb_mean_A.pkl"):
        pb = cc.load_ckpt("pb_mean_A.pkl")
    else:
        pb = cc.pseudobulk(adata, layer=None, agg="mean")
        cc.save_ckpt(pb, "pb_mean_A.pkl")
    n = sum(d.shape[0] for d in pb.values())
    print(f"pseudobulk: {len(pb)} cell types, {n} sample-ct blocks", flush=True)

    if which in ("ebmf", "all") and not cc.ckpt_exists("rep_ebmf_A.pkl"):
        rep = fit_ebmf(pb)
        cc.save_ckpt(rep, "rep_ebmf_A.pkl")
        print(f"EBMF-A rep: {rep.shape}", flush=True)
    if which in ("mofa", "all") and not cc.ckpt_exists("rep_mofa_A.pkl"):
        rep = fit_mofa(pb)
        cc.save_ckpt(rep, "rep_mofa_A.pkl")
        print(f"MOFA-A rep: {rep.shape}", flush=True)
    print(f"regime A done ({time.time()-t0:.0f}s)", flush=True)
