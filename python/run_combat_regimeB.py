"""Regime B fits: MOFAcellulaR author pipeline (TMM + scran HVGs) and B2 ablation
(TMM, all genes). EBMF-B = locked config identical to regime A by construction.

Usage: python run_combat_regimeB.py [mofa|all]
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    adata, st = cc.load_combat()
    print(f"loaded: {adata.shape} ({time.time()-t0:.0f}s)", flush=True)
    if cc.ckpt_exists("pb_sum_B.pkl"):
        pb = cc.load_ckpt("pb_sum_B.pkl")
    else:
        pb = cc.pseudobulk(adata, layer="X_raw_counts", agg="sum")
        cc.save_ckpt(pb, "pb_sum_B.pkl")
    print(f"raw-count pseudobulk: {len(pb)} cell types", flush=True)
    del adata

    from run_combat_regimeA import fit_mofa
    if which in ("mofa", "all") and not cc.ckpt_exists("rep_mofa_B.pkl"):
        rep = fit_mofa(pb, tag="B", hvg=True)
        cc.save_ckpt(rep, "rep_mofa_B.pkl")
        print(f"MOFA-B rep: {rep.shape}", flush=True)
    if which in ("mofa", "all") and not cc.ckpt_exists("rep_mofa_B2.pkl"):
        rep = fit_mofa(pb, tag="B2", hvg=False)
        cc.save_ckpt(rep, "rep_mofa_B2.pkl")
        print(f"MOFA-B2 rep: {rep.shape}", flush=True)
    print(f"regime B done ({time.time()-t0:.0f}s)", flush=True)
