# Extract per-donor per-cell-type cell counts from each cohort h5ad (obs only).
# Output: one CSV per cohort with columns donor_id, cell_type, n_cells.
#
# Usage: python3 composition_counts.py <h5ad_path> <cohort_name>
# Writes CSV to stdout. Reads only the obs group (fast even on large files).
# Drops the "Bad quality" annotation group, consistent with pseudobulking.
import h5py
import numpy as np
import pandas as pd
import sys
import os

h5ad_path, cohort = sys.argv[1], sys.argv[2]

with h5py.File(h5ad_path, "r") as f:
    obs = f["obs"]
    def decode_cat(key):
        g = obs[key]
        if "categories" in g:  # pandas categorical
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                             for c in g["categories"][:]])
            codes = g["codes"][:]
            return cats[codes]
        vals = g[:]
        if vals.dtype.kind in "SU":
            return np.array([c.decode() if isinstance(c, bytes) else str(c) for c in vals])
        return vals.astype(str)

    donor = decode_cat("donor_id")
    ct = decode_cat("annotation")

df = pd.DataFrame({"donor_id": donor, "cell_type": ct})
df = df[df.cell_type != "Bad quality"]
counts = (df.groupby(["donor_id", "cell_type"]).size().reset_index(name="n_cells"))
counts.to_csv(sys.stdout, index=False)
sys.stderr.write(f"{cohort}: {counts.shape[0]} donor-ct pairs, "
                 f"{counts.donor_id.nunique()} donors\n")
