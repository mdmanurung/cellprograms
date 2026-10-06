"""Baseline representations for COMBAT regime A (matched input).

PCA: 30 PCs of the same mean pseudobulk used by EBMF.
Pseudobulk / GroupedPseudobulk: patpy classes on log-normalized genes.
RandomVector: 30-dim gaussian.

Usage: python run_combat_baselines.py
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)
SEED = 67


def main():
    t0 = time.time()
    adata, st = cc.load_combat()
    pb = cc.load_ckpt("pb_mean_A.pkl") if cc.ckpt_exists("pb_mean_A.pkl") else None
    if pb is None:
        pb = cc.pseudobulk(adata, agg="mean")
        cc.save_ckpt(pb, "pb_mean_A.pkl")

    # PCA on the same mean pseudobulk (union of samples, zero-filled blocks)
    if not cc.ckpt_exists("rep_pca_A.pkl"):
        blocks = []
        for ct, df in pb.items():
            z = (df - df.mean(0)) / (df.std(0, ddof=0) + 1e-12)
            z.columns = [f"{ct}::{c}" for c in z.columns]
            blocks.append(z)
        X = pd.concat(blocks, axis=1).fillna(0.0)
        from sklearn.decomposition import PCA
        Z = PCA(n_components=30, random_state=SEED).fit_transform(X.values)
        rep = pd.DataFrame(Z, index=X.index, columns=[f"PC{i+1}" for i in range(30)])
        cc.save_ckpt(rep, "rep_pca_A.pkl")
        print(f"PCA rep: {rep.shape}", flush=True)
    from patpy.tl.sample_representation import Pseudobulk, GroupedPseudobulk, RandomVector
    adata.obs[cc.SAMPLE_KEY] = adata.obs[cc.SAMPLE_KEY].astype(str)

    if not cc.ckpt_exists("rep_pseudobulk_A.pkl"):
        m = Pseudobulk(sample_key=cc.SAMPLE_KEY, cell_group_key=cc.CT_KEY, layer=None, seed=SEED)
        m.prepare_anndata(adata)
        m.calculate_distance_matrix()
        cc.save_ckpt(pd.DataFrame(m.sample_representation, index=m.samples), "rep_pseudobulk_A.pkl")
        print("Pseudobulk done", flush=True)
    if not cc.ckpt_exists("rep_gpb_A.pkl"):
        m = GroupedPseudobulk(sample_key=cc.SAMPLE_KEY, cell_group_key=cc.CT_KEY, layer=None, seed=SEED)
        m.prepare_anndata(adata)
        D = m.calculate_distance_matrix()
        cc.save_ckpt({"distances": D, "samples": list(m.samples)}, "rep_gpb_A.pkl")
        print("GroupedPseudobulk done", flush=True)
    if not cc.ckpt_exists("rep_rand_A.pkl"):
        m = RandomVector(sample_key=cc.SAMPLE_KEY, cell_group_key=cc.CT_KEY, latent_dim=30, seed=SEED)
        m.prepare_anndata(adata)
        m.calculate_distance_matrix()
        cc.save_ckpt(pd.DataFrame(m.sample_representation, index=m.samples), "rep_rand_A.pkl")
        print("RandomVector done", flush=True)
    print(f"baselines done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
