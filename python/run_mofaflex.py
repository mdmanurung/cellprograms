"""MOFA-FLEX comparator (bioFAM/mofaflex, mc-ASTRA engine).

Fits MOFA-FLEX on the same pseudobulk arms as the EBMF tuning search:
views = cell types, samples = donors, shared latent Z (samples x K).
Configs: weight_prior in {Normal (MOFA2-like), Horseshoe (structured
sparsity)}, n_factors=10 (matched to MOFA-B), init_factors='pca', seed=67.

Usage: python run_mofaflex.py <arm_dir> <prefix> <out_csv> [--zdir DIR]
       prefix 'mean_A' loads the pb_mean_A.pkl checkpoint instead.
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import combat_common as cc

warnings.filterwarnings("ignore")

CONFIGS = [{"weight_prior": "Normal", "n_factors": 10},
           {"weight_prior": "Horseshoe", "n_factors": 10}]


def load_arm_or_ckpt(arm_dir, prefix):
    if prefix == "mean_A":
        return cc.load_ckpt("pb_mean_A.pkl")
    from run_tuning_search import load_arm
    return load_arm(arm_dir, prefix)


def fit_mofaflex(pb, cfg, seed=67):
    import torch, anndata as ad
    import mofaflex as mf
    torch.set_num_threads(8)
    views = {ct: ad.AnnData(df.astype(np.float32)) for ct, df in pb.items()}
    mo = mf.ModelOptions(n_factors=int(cfg["n_factors"]),
                         weight_prior=cfg["weight_prior"], init_factors="pca")
    to = mf.TrainingOptions(device="cpu", max_epochs=10000, seed=seed,
                            early_stopper_patience=100,
                            save_path="/workspace/mofaflex_h5")
    os.makedirs("/workspace/mofaflex_h5", exist_ok=True)
    model = mf.MOFAFLEX({"single_group": views}, mo, to)
    Z = model.get_factors()
    Z = Z[list(Z.keys())[0]]
    W = dict(model.get_weights())
    return Z, W, model


if __name__ == "__main__":
    arm_dir, prefix, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    zdir = sys.argv[sys.argv.index("--zdir") + 1] if "--zdir" in sys.argv else "/workspace/mofaflex_z"
    os.makedirs(zdir, exist_ok=True)
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    from run_tuning_search import score_rep
    pb = load_arm_or_ckpt(arm_dir, prefix)
    print(f"arm '{prefix}' loaded: {len(pb)} cell types", flush=True)

    rows, done = [], set()
    if os.path.exists(out_csv):
        d = pd.read_csv(out_csv)
        rows = d.to_dict("records")
        done = set(d.loc[d["knn_disease"].notna(), "config"])

    for i, cfg in enumerate(CONFIGS):
        tag = f"mofaflex_{prefix}_{cfg['weight_prior']}_k{cfg['n_factors']}"
        if tag in done:
            print(f"skip {tag} (done)", flush=True)
            continue
        t0 = time.time()
        try:
            Z, W, model = fit_mofaflex(pb, cfg)
            Z.index = Z.index.astype(str)
            s = score_rep(Z, st)
            rt = round(time.time() - t0)
            rows.append({**cfg, **s, "runtime_s": rt, "config": tag})
            Z.to_csv(f"{zdir}/{tag}_Z.csv")
            pd.DataFrame(W.items()).to_csv(f"{zdir}/{tag}_Wshapes.csv", index=False)
            pd.DataFrame(rows).to_csv(out_csv, index=False)
            print(f"{tag}: disease={s['knn_disease']:.3f} who={s['knn_who']:.3f} "
                  f"K={s['n_factors']} ({rt}s)", flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            rows.append({**cfg, "error": str(e)[:200], "config": tag})
            pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("MOFAFLEX DONE", flush=True)
