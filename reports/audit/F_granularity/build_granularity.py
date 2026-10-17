"""Workstream F (PILOT): build three input granularities from one cohort + cell type.

Granularities (all outputs are RAW COUNTS, genes restricted to top-2000 HVG):
  1. pseudobulk  : donors x genes, counts summed per donor within the cell type
  2. metacell    : metacells x genes, k-means within donor on PCs, gamma ~= 25 cells/metacell,
                   counts summed per metacell (pilot simplification; MR-3 will use SEACells + mcRigor)
  3. singlecell  : <=30k cells stratified by donor (proportional), identical subsample for all methods

Usage:
  python build_granularity.py inspect --h5ad data/UCSF_reannotated.h5ad
  python build_granularity.py build   --h5ad data/UCSF_reannotated.h5ad --celltype "..." --outdir gran
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

SEED = 20261017
GAMMA = 25
N_HVG = 2000
MAX_CELLS = 30000
N_PCS = 30

DONOR_CANDS = ["donor", "donor_id", "sample", "sample_id", "subject", "subject_id", "patient", "patient_id"]
TYPE_CANDS = ["author_celltype", "celltype", "cell_type", "CellType", "celltype_major", "annotation", "cell_ontology_class", "manual_annotation"]
AGE_CANDS = ["age", "age_years", "age_yrs", "Age"]


def find_col(cols, cands):
    for c in cands:
        if c in cols:
            return c
    for c in cols:
        lc = str(c).lower()
        for cand in cands:
            if cand.lower() in lc:
                return c
    return None


def counts_of(adata):
    if adata.raw is not None and adata.raw.X is not None:
        X = adata.raw.X
        genes = pd.Index(adata.raw.var_names.astype(str))
    else:
        X = adata.X
        genes = pd.Index(adata.var_names.astype(str))
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    return X.tocsc(), genes


def inspect(h5ad):
    a = ad.read_h5ad(h5ad, backed="r")
    print("shape:", a.shape)
    print("obs columns:", list(a.obs.columns))
    print("raw is not None:", a.raw is not None)
    if a.raw is not None:
        print("raw shape:", a.raw.shape)
    a.file.close()


def build(h5ad, celltype, outdir, type_col=None, donor_col=None, age_col=None):
    os.makedirs(outdir, exist_ok=True)
    a = ad.read_h5ad(h5ad)
    X, genes = counts_of(a)
    obs = a.obs.copy()
    obs.index = obs.index.astype(str)

    type_col = type_col or find_col(obs.columns, TYPE_CANDS)
    donor_col = donor_col or find_col(obs.columns, DONOR_CANDS)
    age_col = age_col or find_col(obs.columns, AGE_CANDS)
    assert type_col is not None and donor_col is not None, "could not find cell type / donor columns"
    print(f"using type_col={type_col!r} donor_col={donor_col!r} age_col={age_col!r}")

    types = obs[type_col].astype(str)
    if celltype is None:
        print(types.value_counts().to_string())
        sys.exit(0)
    mask = (types == celltype).values
    # QC: drop cells flagged bad quality in minor_celltypes if present
    if "minor_celltypes" in obs.columns:
        bad = (obs["minor_celltypes"].astype(str) == "Bad quality").values
        print(f"QC: dropping {int((mask & bad).sum())} Bad-quality cells within type")
        mask = mask & ~bad
    print(f"cell type {celltype!r}: {mask.sum()} cells")
    assert mask.sum() > 0, f"no cells found for type {celltype!r} in column {type_col!r}; levels: {sorted(types.unique())[:20]}"
    X = X[mask, :]
    obs = obs.loc[mask].copy()

    # drop zero-count cells
    lib = np.asarray(X.sum(axis=1)).ravel()
    keep = lib > 0
    if (~keep).sum() > 0:
        print(f"dropping {(~keep).sum()} cells with 0 counts")
    X = X[keep, :]
    obs = obs.iloc[np.where(keep)[0]]
    donor_labels = obs[donor_col].astype(str).values
    uniq_donors = pd.unique(donor_labels)
    print(f"donors: {len(uniq_donors)}")

    # HVG on log1p-CP10K of this cell type
    ad_sub = ad.AnnData(X.tocsr(), obs=obs[[donor_col]].copy(), var=pd.DataFrame(index=genes))
    sc.pp.normalize_total(ad_sub, target_sum=1e4)
    sc.pp.log1p(ad_sub)
    try:
        sc.pp.highly_variable_genes(ad_sub, n_top_genes=N_HVG, flavor="seurat_v3", subset=False)
    except Exception as e:
        print(f"seurat_v3 HVG failed ({e}); falling back to cell_ranger on log data")
        sc.pp.highly_variable_genes(ad_sub, n_top_genes=N_HVG, flavor="cell_ranger", subset=False)
    hv = ad_sub.var["highly_variable"].values
    hvg = list(genes[hv])
    print(f"HVGs: {len(hvg)}")
    gidx = pd.Index(genes).get_indexer(hvg)
    Xh = X[:, gidx].tocsc()  # counts, cells x hvg

    # PCA space for metacell clustering (log-normalized, HVG only)
    L = ad_sub[:, hv].X
    L = np.asarray(L.todense()) if sp.issparse(L) else np.asarray(L)
    n_pcs = int(min(N_PCS, L.shape[0] - 1, L.shape[1] - 1))
    Z = PCA(n_components=n_pcs, random_state=SEED).fit_transform(L)

    donor_index = {}
    for i, d in enumerate(donor_labels):
        donor_index.setdefault(d, []).append(i)

    # ---- granularity 1: pseudobulk (donors x genes) ----
    pb = np.vstack([np.asarray(Xh[donor_index[d]].sum(axis=0)).ravel() for d in uniq_donors])
    print("pseudobulk:", pb.shape)

    # ---- granularity 2: metacells within donor x cell type ----
    mc_rows, mc_donor, mc_size = [], [], []
    for d in uniq_donors:
        idx = np.asarray(donor_index[d])
        k = max(1, int(round(len(idx) / GAMMA)))
        if k == 1:
            labs = np.zeros(len(idx), dtype=int)
        else:
            labs = KMeans(n_clusters=k, n_init=10, random_state=SEED).fit_predict(Z[idx])
        for lab in np.unique(labs):
            members = idx[labs == lab]
            mc_rows.append(np.asarray(Xh[members].sum(axis=0)).ravel())
            mc_donor.append(d)
            mc_size.append(int(len(members)))
    M = np.vstack(mc_rows)
    print(f"metacells: {M.shape} (mean gamma {np.mean(mc_size):.1f})")

    # ---- granularity 3: single-cell stratified subsample ----
    rng = np.random.default_rng(SEED)
    n_total = len(donor_labels)
    take = []
    for d in uniq_donors:
        idx = np.asarray(donor_index[d])
        n_take = min(len(idx), max(1, int(round(len(idx) * MAX_CELLS / n_total))))
        take.extend(rng.choice(idx, size=n_take, replace=False))
    take = np.sort(np.asarray(take))
    if len(take) > MAX_CELLS:
        take = np.sort(rng.choice(take, size=MAX_CELLS, replace=False))
    S = Xh[take].toarray()
    sc_donor = donor_labels[take]
    print("single-cell:", S.shape)

    # ---- write everything ----
    from scipy.io import mmwrite
    spec = {
        "pseudobulk": (pb, np.array([f"donor_{d}" for d in uniq_donors])),
        "metacell": (M, np.array([f"mc{i}" for i in range(M.shape[0])])),
        "singlecell": (S, np.array([f"cell{i}" for i in range(S.shape[0])])),
    }
    meta = {
        "pseudobulk": pd.DataFrame({"id": spec["pseudobulk"][1], "donor": uniq_donors}),
        "metacell": pd.DataFrame({"id": spec["metacell"][1], "donor": mc_donor, "n_cells": mc_size}),
        "singlecell": pd.DataFrame({"id": spec["singlecell"][1], "donor": sc_donor}),
    }
    for name, (mat, ids) in spec.items():
        d = os.path.join(outdir, name)
        os.makedirs(d, exist_ok=True)
        mmwrite(os.path.join(d, "counts.mtx"), sp.csc_matrix(mat))
        meta[name].to_csv(os.path.join(d, "obs_meta.tsv"), sep="\t", index=False)
        with open(os.path.join(d, "genes.tsv"), "w") as f:
            f.write("\n".join(hvg))
    recipe = dict(cohort=os.path.basename(h5ad), celltype=celltype, type_col=type_col, donor_col=donor_col,
                  age_col=age_col, seed=SEED, gamma=GAMMA, n_hvg=N_HVG, max_cells=MAX_CELLS, n_pcs=n_pcs,
                  n_donors=int(pb.shape[0]), n_metacells=int(M.shape[0]), n_singlecell=int(S.shape[0]),
                  n_cells_type=int(n_total), hvg=hvg)
    with open(os.path.join(outdir, "recipe.json"), "w") as f:
        json.dump(recipe, f, indent=1)
    if age_col is not None:
        obs.groupby(obs[donor_col].astype(str))[age_col].first().to_csv(
            os.path.join(outdir, "donor_age.tsv"), sep="\t")
    print("done")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["inspect", "build"])
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--celltype", default=None)
    ap.add_argument("--outdir", default="gran")
    ap.add_argument("--type-col", default=None)
    ap.add_argument("--donor-col", default=None)
    ap.add_argument("--age-col", default=None)
    args = ap.parse_args()
    if args.mode == "inspect":
        inspect(args.h5ad)
    else:
        build(args.h5ad, args.celltype, args.outdir, args.type_col, args.donor_col, args.age_col)
