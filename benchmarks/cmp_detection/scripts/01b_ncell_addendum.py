#!/usr/bin/env python3
## 01b_ncell_addendum.py <h5ad> <annotation_csv> <semi_dir>
## Writes <semi_dir>/<type>/ncell.csv (donor,ncell) from obs + annotation only
## (no X read; ~1 min). gen_semi() consumes b$cts[[ct]]$ncell.
import sys, os
import h5py
import numpy as np
import pandas as pd

# MAJOR mapping + _h5_col inlined (the 01 script defines them inside main())
MAJOR = {
    "B": ["B_naive", "B_non-switched_memory", "B_switched_memory", "B_immature", "B_exhausted"],
    "Plasma": ["Plasma_cell_IgA", "Plasma_cell_IgG", "Plasma_cell_IgM", "Plasmablast"],
    "CD4_T": ["CD4.Naive", "CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Tfh", "CD4.Th1", "CD4.Th17", "CD4.Th2", "CD4.Prolif", "Treg"],
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

def _h5_col(group, key):
    g = group[key]
    if isinstance(g, h5py.Group) and "categories" in g and "codes" in g:
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"]])
        codes = np.asarray(g["codes"])
        return np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "NA").astype(str)
    d = np.asarray(g[()])
    if d.dtype.kind in "iu" and "__categories" in group and key in group["__categories"]:
        cats = np.array([c.decode() if isinstance(c, bytes) else c
                         for c in group["__categories"][key][()]])
        return np.where(d >= 0, cats[np.clip(d, 0, len(cats) - 1)], "NA").astype(str)
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in d.reshape(-1)]) \
        if d.dtype.kind in "OSU" else d.reshape(-1)

h5ad_path, annot_path, semi_dir = sys.argv[1], sys.argv[2], sys.argv[3]
f = h5py.File(h5ad_path, "r")
obs_g = f["obs"]
idx_key = next((k for k in ("_index", "index", "covid_index") if k in obs_g), None)
obs_ids = _h5_col(obs_g, idx_key)
donors = _h5_col(obs_g, "patient_id").astype(str)
ann = pd.read_csv(annot_path, index_col=0)
fine2major = {fl: mj for mj, fines in MAJOR.items() for fl in fines}
majors = np.array([fine2major.get(b, "EXCLUDE") for b in ann["full_clustering"]])
# align annotation rows to obs ids
pos = pd.Series(np.arange(len(obs_ids)), index=obs_ids)
ann_bc = pd.Index(ann.index.astype(str))
sel = ann_bc.isin(pos.index).values
majors = majors[sel]
cell_pos = pos.loc[ann_bc[sel]].values
donor_sel = donors[cell_pos]
for mj in sorted(set(MAJOR)):
    m = majors == mj
    if not m.any():
        continue
    d = donor_sel[m]
    vc = pd.Series(d).value_counts().sort_index()
    out = os.path.join(semi_dir, mj)
    os.makedirs(out, exist_ok=True)
    pd.DataFrame({"donor": vc.index, "ncell": vc.values}).to_csv(
        os.path.join(out, "ncell.csv"), index=False)
    print(f"{mj}: {len(vc)} donors, {int(vc.sum())} cells")
