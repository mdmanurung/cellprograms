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
    "HSC": ["HSC"],
    "RBC": ["RBC"],
}
EXCLUDED = ["B_malignant"]  # malignant cells are not a normal donor program

MIN_CELLS = 50


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

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    print("obs columns:", list(obs.columns))
    # donor/site/sex/status columns -- verify names on first run
    # columns PINNED from the actual E-MTAB-10026 h5ad obs (verified 2026-10-22);
    # heuristics kept only as fallback
    donor_col = next((c for c in obs.columns if c.lower() in ("patient_id", "donor_id", "donor")), None) or \
        next(c for c in obs.columns if "sample" in c.lower())
    site_col = next(c for c in obs.columns if "site" in c.lower() or c.lower() in ("institute", "hospital"))
    sex_col = next((c for c in obs.columns if c.lower() == "sex"), None)
    status_col = next((c for c in obs.columns if c.lower() == "status_on_day_collection_summary"), None) or \
        next(c for c in obs.columns if "status" in c.lower())

    bc2major = ann["major"].to_dict()
    keep_mask = np.array([b in bc2major for b in obs_names(adata)])
    print(f"cells in annotation: {keep_mask.sum()} / {len(keep_mask)}")

    # counts: prefer .raw (raw counts) if X is not integer-valued
    X = adata.X if adata.raw is None else adata.raw.X
    counts_are_int = is_integerish(X[:200] if sparse.issparse(X) else X[:200])
    if adata.raw is not None and not counts_are_int:
        print("X is not counts; using .raw.X")
    genes = list(adata.raw.var_names if (adata.raw is not None and not counts_are_int) else adata.var_names)

    majors = np.array([bc2major.get(b, "EXCLUDE") for b in obs_names(adata)])
    donors_all = np.array(obs[donor_col].astype(str))
    site = np.array(obs[site_col].astype(str))
    sex = np.array(obs[sex_col].astype(str)) if sex_col else ["NA"] * len(obs)
    status = np.array(obs[status_col].astype(str))

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
        with gzip.open(os.path.join(tdir, "Y.mtx.gz"), "wb") as fh:
            scio.mmwrite(fh, Y)
        with open(os.path.join(tdir, "donors.txt"), "w") as fh:
            fh.write("\n".join(ud) + "\n")
        pd.DataFrame({"gene": genes, "cpm_mean": cpm_mean}).to_csv(
            os.path.join(tdir, "cpm_mean.csv"), index=False)
        report.append(f"{major}: {len(ud)} donors, {len(sel)} cells, {Y.shape[1]} genes")
        dm = pd.DataFrame({
            "donor_id": ud,
            "site": [np.unique(site[donors_all[sel] == d])[0] for d in ud],
            "sex": [np.unique(sex[donors_all[sel] == d])[0] for d in ud],
            "Status_on_day_collection_summary": [
                np.unique(status[donors_all[sel] == d])[0] for d in ud],
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
