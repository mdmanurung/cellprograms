"""COMBAT benchmark: shared data loading, targets, pseudobulk, checkpoints.

Regime A: matched input — log-normalized pseudobulk, all 3000 genes, same
samples/cell types for every method.
Regime B: author-recommended — raw-count pseudobulk -> TMM (scale 1e6) ->
per-view scran HVGs (B) or all genes (B2) -> MOFAcellulaR/MOFA2.
"""
import os, sys, time, json, warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

WORK = "/workspace/combat_benchmark"
CKPT = "/mnt/shared-workspace/combat_benchmark"
os.makedirs(WORK, exist_ok=True)
os.makedirs(CKPT, exist_ok=True)

SAMPLE_KEY = "scRNASeq_sample_ID"
CT_KEY = "Annotation_major_subset"
MIN_CELLS_PER_SAMPLE_CT = 10  # matched-input filter: drop (sample, ct) pairs below this


def load_combat():
    """Load COMBAT, apply exclusions, return (adata, sample_table)."""
    adata = ad.read_h5ad("/workspace/combat_processed.h5ad")
    meta = ad.read_h5ad("/workspace/combat_meta.h5ad").obs.copy()
    meta.index = meta.index.astype(str)

    obs = adata.obs.copy()
    obs[SAMPLE_KEY] = obs[SAMPLE_KEY].astype(str)
    # drop QC-flagged samples
    excl = set(meta.index[meta["potential_exclude"].astype(bool)])
    keep_samples = set(meta.index)  # also drops the 1 sample absent from meta
    mask = obs[SAMPLE_KEY].isin(keep_samples - excl)
    adata = adata[mask.values]
    obs = adata.obs.copy()

    # sample-level table
    st = meta.loc[sorted(obs[SAMPLE_KEY].unique())].copy()
    dc = st["DiseaseClassification"].astype(str).str.split(";").str[0]
    st["disease"] = dc.replace({"NA": "Healthy"})
    st["source"] = st["Source"].astype(str)
    st["who_ordinal"] = st["WHO_ordinal_at_sample"]
    st["outcome"] = st["Outcome"]  # reverse-coded vs WHO (6=best, 1=worst)
    st["institute"] = st["Institute"].astype(str)
    st["pool"] = st["Pool_ID"].astype(str)
    return adata, st


def pseudobulk(adata, layer=None, agg="sum", ct_key=CT_KEY, sample_key=SAMPLE_KEY):
    """Samples x (celltype x genes) pseudobulk. Returns dict ct -> DataFrame.

    agg="sum": raw-count sums (regime B, TMM input). agg="mean": per-sample
    means (regime A, matched with the adapter's pseudobulk).
    """
    X = adata.layers[layer] if layer else adata.X
    if sp.issparse(X):
        X = X.tocsr()
    obs = adata.obs
    genes = adata.var_names.astype(str).tolist()
    out = {}
    for ct in obs[ct_key].cat.categories if hasattr(obs[ct_key], "cat") else obs[ct_key].unique():
        m = (obs[ct_key] == ct).values
        if m.sum() == 0:
            continue
        sub_obs = obs.loc[m]
        sub_X = X[m]
        # indicator matrix: cells x samples
        samples = sub_obs[sample_key].astype(str).values
        uniq = pd.Index(sorted(set(samples)))
        ind = sp.csr_matrix(
            (np.ones(len(samples)), (np.arange(len(samples)), uniq.get_indexer(samples))),
            shape=(len(samples), len(uniq)),
        )
        pb = ind.T @ sub_X  # sums per sample
        counts = np.asarray(ind.sum(axis=0)).ravel()
        pb = pb.toarray() if sp.issparse(pb) else pb
        if agg == "mean":
            pb = pb / np.maximum(counts, 1)[:, None]
        df = pd.DataFrame(pb, index=uniq.astype(str), columns=genes)
        # matched-input filter: min cells per (sample, ct)
        keep = counts >= MIN_CELLS_PER_SAMPLE_CT
        df = df.loc[uniq.astype(str)[keep]]
        if len(df) > 0:
            out[str(ct)] = df
    return out


def save_ckpt(obj, name):
    import pickle
    with open(f"{CKPT}/{name}", "wb") as f:
        pickle.dump(obj, f)
    # verify
    assert os.path.getsize(f"{CKPT}/{name}") > 0
    print(f"checkpointed {name} ({os.path.getsize(f'{CKPT}/{name}')/1e6:.1f} MB)", flush=True)


def load_ckpt(name):
    import pickle
    with open(f"{CKPT}/{name}", "rb") as f:
        return pickle.load(f)


def ckpt_exists(name):
    return os.path.exists(f"{CKPT}/{name}")
