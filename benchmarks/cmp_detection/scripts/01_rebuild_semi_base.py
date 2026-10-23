#!/usr/bin/env python3
"""01_rebuild_semi_base.py -- rebuild the Stephenson pseudobulk substrate.

Registered in PREREG_mr2-cnmf-v1.md (deviation 1): donor x gene COUNT
pseudobulk per cell type from ArrayExpress E-MTAB-10026
(covid_portal_210320_with_raw.h5ad + annotation_200112.csv).

Frozen 12-major mapping (below). Donors with < 50 cells of a type are
dropped FOR THAT TYPE. Outputs (SEMI_DIR, default
benchmarks/cmp_detection/data/stephenson_pb/):
  <type>/Y.mtx.gz        donor x gene counts (MatrixMarket, donor-major)
  <type>/donors.txt      donor ids (row order of Y)
  <type>/cpm_mean.csv    donor x gene CPM (size-factor 1e6)
  donor_meta.csv         donor_id, site, sex, Status_on_day_collection_summary
  rebuild_report.txt     per-type donor/cell counts + mapping coverage assert

Usage: python3 01_rebuild_semi_base.py <h5ad> <annotation_csv> [semi_dir]
"""
import gzip
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as scio
from scipy import sparse

# ---- frozen 12-major mapping (PREREG_mr2-cnmf-v1.md deviation 1) ------------
MAJOR = {
    "B": ["B_naive", "B_non-switched_memory", "B_switched_memory",
          "B_immature", "B_exhausted"],
    "Plasma": ["Plasma_cell_IgA", "Plasma_cell_IgG", "Plasma_cell_IgM",
               "Plasmablast"],
    "CD4_T": ["CD4.Naive", "CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Tfh",
              "CD4.Th1", "CD4.Th17", "CD4.Th2", "CD4.Prolif", "Treg"],
    "CD8_T": ["CD8.Naive", "CD8.EM", "CD8.TE", "CD8.Prolif"],
    "NK": ["NK_16hi", "NK_56hi", "NK_prolif"],
    "gdT": ["gdT", "MAIT", "NKT", "ILC1_3", "ILC2"],
    "CD14_mono": ["CD14_mono", "CD83_CD14_mono", "Mono_prolif"],
    "CD16_mono": ["CD16_mono", "C1_CD16_mono"],
    "DC": ["DC1", "DC2", "DC3", "pDC", "DC_prolif", "ASDC"],
    "Platelets": ["Platelets"],
    "HSC": ["HSC", "HSC_CD38pos", "HSC_myeloid", "HSC_CD38neg",
            "HSC_prolif", "HSC_erythroid", "HSC_MK"],
    "RBC": ["RBC"],
}
EXCLUDED = ["B_malignant"]  # malignant cells are not a normal donor program

MIN_CELLS = 50
# Cell types with fewer than MIN_DONORS donors cannot support the registered
# CMP-detection K (max_factors = 30 needs n >= K + 2 = 32 under the C4
# residual-df guard, R/ebmf.R). The lost original semi_base predates that
# guard; the rebuild drops such types so every retained type is K30-fittable.
MIN_DONORS = 32


def main(h5ad_path, annot_path, semi_dir):
    ann = pd.read_csv(annot_path, index_col=0)
    fine2major = {}
    for major, fines in MAJOR.items():
        for f in fines:
            fine2major[f] = major
    unknown = set(ann["full_clustering"]) - set(fine2major) - set(EXCLUDED)
    assert not unknown, f"mapping does not cover fine types: {unknown}"
    ann["major"] = [fine2major.get(f, "EXCLUDE") for f in ann["full_clustering"]]
    ann = ann[ann["major"] != "EXCLUDE"]

    # ---- direct h5py read (the as-run approach) --------------------------------
    # anndata's eager .uns read crashes on this OLD-format h5ad (encoding-version
    # 0.1.0: no raw, no shape dataset in X, categoricals as int codes +
    # obs/__categories, and a corrupt gzip block in X/indices). Read via h5py.
    import h5py
    f = h5py.File(h5ad_path, "r")
    obs_g = f["obs"]

    def h5col(group, key):
        # old-format categorical: int codes + labels in obs/__categories/<col>
        d = np.asarray(group[key][()])
        if d.dtype.kind in "iu" and "__categories" in group and key in group["__categories"]:
            cats = np.array([c.decode() if isinstance(c, bytes) else c
                             for c in group["__categories"][key][()]])
            return np.where(d >= 0, cats[np.clip(d, 0, len(cats) - 1)], "NA").astype(str)
        if d.dtype.kind in "OSU":
            return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in d.reshape(-1)])
        return d.reshape(-1)

    # obs index: this file stores barcodes in `covid_index` (verified 2026-10-22);
    # fall back through the usual keys
    idx_key = next((k for k in ("covid_index", "_index", "index") if k in obs_g), None)
    assert idx_key, f"no obs index key; have {list(obs_g.keys())}"
    barcodes = h5col(obs_g, idx_key)
    # columns PINNED from the actual E-MTAB-10026 h5ad obs (verified 2026-10-22)
    donors_all = h5col(obs_g, "patient_id").astype(str)
    site = h5col(obs_g, "Site").astype(str)
    sex = h5col(obs_g, "Sex").astype(str)
    status = h5col(obs_g, "Status_on_day_collection_summary").astype(str)

    bc2major = ann["major"].to_dict()
    majors = np.array([bc2major.get(b, "EXCLUDE") for b in barcodes])
    keep_mask = majors != "EXCLUDE"
    print(f"cells in annotation: {keep_mask.sum()} / {len(keep_mask)}")

    # X: float32 CSR, no shape dataset (old format) -> derive from obs/var sizes
    Xg = f["X"] if isinstance(f["X"], h5py.Group) else None
    var_g = f["var"]
    vidx = next((k for k in ("_index", "index") if k in var_g), None)
    assert vidx, f"no var index key; have {list(var_g.keys())}"
    genes = [g.decode() if isinstance(g, bytes) else str(g) for g in var_g[vidx][()]]
    n_cells, n_genes = len(barcodes), len(genes)
    print(f"X: {n_cells} x {n_genes}")

    CHUNK = 10_000_000  # large reads trip the gzip filter; chunk at 10M elements
    nnz = Xg["data"].shape[0]
    # int32 indices: n_genes < 2^31 and this halves peak memory (7.7 GB -> 3.9 GB)
    data = np.zeros(nnz, dtype=np.float32)
    indices = np.zeros(nnz, dtype=np.int32)
    for s in range(0, nnz, CHUNK):
        e = min(s + CHUNK, nnz)
        data[s:e] = Xg["data"][s:e]
        try:
            indices[s:e] = Xg["indices"][s:e]
        except OSError:
            # corrupt gzip block (~indices[559,150,000:559,250,000)): zero these
            # CSR entries (data=0 contributes nothing to donor sums); disclosed
            # deviation in PREREG_mr2-cnmf-v1.md
            try:
                idx_chunk = Xg["indices"][s:e]
                idx_chunk = np.nan_to_num(np.mod(idx_chunk, max(n_genes, 1))).astype(np.int32)
            except OSError:
                idx_chunk = np.zeros(e - s, dtype=np.int32)
            data[s:e] = 0.0
            indices[s:e] = idx_chunk
            print(f"WARNING: zeroed corrupt indices block [{s},{e})")
    indptr = np.asarray(Xg["indptr"][()], dtype=np.int64)
    X = sparse.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))
    f.close()
    counts_are_int = is_integerish(X[:200] if sparse.issparse(X) else X[:200])
    print("X integer-valued (counts):", bool(counts_are_int))

    os.makedirs(semi_dir, exist_ok=True)
    report = []
    donor_meta = None
    for major in sorted(set(MAJOR)):
        sel = np.where(majors == major)[0]
        if not len(sel):
            continue
        # <50 cells/donor dropped for this type
        sub_donors = donors_all[sel]
        vc = pd.Series(sub_donors).value_counts()
        ok_donors = vc[vc >= MIN_CELLS].index
        if len(ok_donors) < MIN_DONORS:
            report.append(f"{major}: dropped, only {len(ok_donors)} donors with >= {MIN_CELLS} cells (< MIN_DONORS={MIN_DONORS}; cannot support K=30)")
            continue
        sel = sel[np.isin(sub_donors, ok_donors)]
        if not len(sel):
            report.append(f"{major}: 0 donors after >= {MIN_CELLS} cell filter")
            continue
        Xt = X[sel]
        if sparse.issparse(Xt):
            Xt = Xt.tocsr()
        ud = np.unique(donors_all[sel])
        D = sparse.csr_matrix(
            (np.ones(len(sel)), (np.searchsorted(ud, donors_all[sel]), np.arange(len(sel)))),
            shape=(len(ud), len(sel)))
        Y = (D @ Xt).tocsr()          # donor x gene counts
        lib = np.asarray(Y.sum(1)).ravel()
        cpm = Y.multiply((1e6 / np.maximum(lib, 1))[:, None]).tocsr()
        cpm_mean = np.asarray(cpm.mean(0)).ravel()  # type-level mean CPM profile
        tdir = os.path.join(semi_dir, major)
        os.makedirs(tdir, exist_ok=True)
        # Y.mtx.gz stores log2(CPM + 1): the sim_cmp.R contract is "Y logCPM"
        # (the ambient transform inverts it via 2^Y - 1, and the planting
        # amplitude msd is a log-scale sd). The cNMF arm recovers CPM-scale
        # values via 2^Y - 1 (Addendum A as amended).
        Ylog = np.log2(cpm.toarray() + 1) if sparse.issparse(cpm) else np.log2(cpm + 1)
        with gzip.open(os.path.join(tdir, "Y.mtx.gz"), "wb") as fh:
            scio.mmwrite(fh, sparse.csr_matrix(Ylog))
        with open(os.path.join(tdir, "donors.txt"), "w") as fh:
            fh.write("\n".join(ud) + "\n")
        pd.DataFrame({"gene": genes, "cpm_mean": cpm_mean}).to_csv(
            os.path.join(tdir, "cpm_mean.csv"), index=False)
        report.append(f"{major}: {len(ud)} donors, {len(sel)} cells, {Y.shape[1]} genes")
        dm = pd.DataFrame({
            "donor_id": ud,
            "site": [np.unique(site[sel][donors_all[sel] == d])[0] for d in ud],
            "sex": [np.unique(sex[sel][donors_all[sel] == d])[0] for d in ud],
            "Status_on_day_collection_summary": [
                np.unique(status[sel][donors_all[sel] == d])[0] for d in ud],
        })
        dm = dm.rename(columns={"donor_id": "donor_id"})
        if donor_meta is None:
            donor_meta = dm
        else:
            donor_meta = donor_meta.merge(dm, on="donor_id", how="outer", suffixes=("", "_dup"))
            for c in ("site", "sex", "Status_on_day_collection_summary"):
                if c + "_dup" in donor_meta:
                    donor_meta[c] = donor_meta[c].fillna(donor_meta[c + "_dup"])
                    donor_meta = donor_meta.drop(columns=[c + "_dup"])
    donor_meta = donor_meta.sort_values("donor_id")
    donor_meta.to_csv(os.path.join(semi_dir, "donor_meta.csv"), index=False)
    with open(os.path.join(semi_dir, "rebuild_report.txt"), "w") as fh:
        fh.write("\n".join(report) + "\n")
    print("\n".join(report))


def obs_names(adata):
    return np.array(adata.obs_names.astype(str))


def is_integerish(x):
    v = np.asarray(x.todense() if sparse.issparse(x) else x).ravel()
    v = v[v > 0]
    return bool(len(v)) and np.allclose(v, np.round(v))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else
         os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "stephenson_pb"))
