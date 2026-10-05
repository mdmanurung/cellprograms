"""Per-cell-type EBMF fit isolation: which cell type breaks the ebnm solver?"""
import os

os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")

import anndata as ad
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from scipy import sparse as sp

adata = ad.read_h5ad("/workspace/stephenson_processed.h5ad")
adata = adata[(adata.obs["Status"] != "Non_covid").values]
X = adata.X
obs = pd.DataFrame({"sample": adata.obs["sample_id"].astype(str).values,
                    "cell_type": adata.obs["author_cell_type"].astype(str).values},
                   index=adata.obs_names.astype(str))
genes = np.asarray(adata.var_names.astype(str))
sample_labels = sorted(obs["sample"].unique())
sample_pos = {s: i for i, s in enumerate(sample_labels)}

pseudobulk = {}
for ct in sorted(obs["cell_type"].unique()):
    mask = (obs["cell_type"] == ct).values
    rows = np.array([sample_pos[s] for s in obs["sample"].values[mask]])
    n_s = len(sample_labels)
    A = sp.csr_matrix((np.ones(mask.sum()), (rows, np.arange(mask.sum()))), shape=(n_s, mask.sum()))
    counts = np.bincount(rows, minlength=n_s)
    agg = (A @ X[mask]).toarray() / np.maximum(counts, 1)[:, None]
    keep = counts > 0
    pseudobulk[ct] = pd.DataFrame(agg[keep], index=np.array(sample_labels)[keep], columns=genes)

ro.r('.libPaths(c("/workspace/.Rlib", .libPaths())); suppressMessages(library(cellprograms))')
from rpy2.robjects import numpy2ri, pandas2ri

with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
    for ct, mat in pseudobulk.items():
        ro.globalenv["Y"] = mat.values
        ro.globalenv["obsid"] = ro.StrVector(mat.index.tolist())
        ro.globalenv["gn"] = ro.StrVector(mat.columns.tolist())
        rcode = (
            "dimnames(Y) <- list(obsid, gn)\n"
            "res <- tryCatch({\n"
            "  cpd <- as_cell_program_data(list(CT = Y))\n"
            "  fit <- fit_celltype_programs(cpd, seed = 67)\n"
            '  paste0("OK_K=", ncol(fit$scores$CT))\n'
            "}, error = function(e) conditionMessage(e))\n"
            "res\n"
        )
        ok = str(ro.r(rcode)[0])
        print(ct, "->", ok[:100], flush=True)
print("DEBUG DONE")
