#!/usr/bin/env python3
"""Phase 1 cohort audit: inspect a portal *_reannotated.h5ad without loading expression.

Usage: python3 audit_h5ad.py <file.h5ad> [--obs-sample N]

Prints a structured audit: matrix encodings, counts-vs-normalized check for .X
and .raw/X, gene ID format, .obs columns with sample values for identity/
demographics columns, cell-type annotation candidates, donor counts.
"""
import sys
import json
import h5py
import numpy as np

IDENTITY_HINTS = ["donor", "sample", "patient", "subject", "participant", "demux", "barcode"]
DEMO_HINTS = ["age", "sex", "gender", "cmv", "ancestry", "ethnic", "race", "cohort", "dataset", "study", "donor"]
CT_HINTS = ["cell_type", "celltype", "harmonized", "annotation", "compartment", "lineage", "broad", "fine"]


def enc_info(node):
    if node is None:
        return None
    if isinstance(node, h5py.Group):
        out = {"storage": "group keys=" + str(sorted(node.keys()))}
        if "X" in node:
            out["X"] = enc_info(node["X"])
        return out
    out = {"shape": tuple(int(s) for s in node.shape), "dtype": str(node.dtype)}
    if isinstance(node, h5py.Dataset):
        out["storage"] = "dense"
    else:
        out["storage"] = "sparse:" + node.attrs.get("encoding-type", "?")
        for k in ("data", "indices", "indptr"):
            if k in node:
                out[f"{k}_dtype"] = str(node[k].dtype)
    return out


def slice_values(node, n=200_000):
    """Read a small sample of values from a sparse/csc dataset group or dense ds."""
    if node is None:
        return None
    if isinstance(node, h5py.Dataset):
        flat = node[:].ravel()
        return flat[flat != 0][:n]
    data = node["data"]
    take = min(n, data.shape[0])
    return data[:take]


def audit(path, obs_sample=8):
    with h5py.File(path, "r") as f:
        print(f"== {path}")
        print("top-level:", sorted(f.keys()))
        for key in ("X", "raw", "raw/X", "obsp", "varp", "layers", "uns"):
            node = f.get(key)
            if key in ("obsp", "varp", "layers") and node is not None:
                print(f"  {key}: {sorted(node.keys())}")
            elif node is not None:
                print(f"  {key}: {enc_info(node)}")

        # counts vs normalized for X and raw/X
        for key in ("X", "raw/X"):
            vals = slice_values(f.get(key))
            if vals is None or len(vals) == 0:
                print(f"  [{key}] no values sampled")
                continue
            v = np.asarray(vals, dtype=np.float64)
            int_like = np.allclose(v, np.round(v))
            print(f"  [{key}] nonzero-sample: n={len(v)} min={v.min():.3f} max={v.max():.3f} "
                  f"mean={v.mean():.3f} integer_like={bool(int_like)}")

        # obs
        obs = f["obs"]
        cols = list(obs.keys())
        n_cells = f["obs"].get("_index")
        n_cells = n_cells.shape[0] if n_cells is not None else None
        print(f"  obs columns ({len(cols)}): {cols}")
        for c in cols:
            g = obs[c]
            if isinstance(g, h5py.Group) and "categories" in g:  # categorical
                cats = g["categories"][:]
                cats = [x.decode() if isinstance(x, bytes) else str(x) for x in cats[:40]]
                print(f"    {c}: categorical n={len(cats)} e.g. {cats[:12]}")
            else:
                d = g if isinstance(g, h5py.Dataset) else g.get("data")
                if d is not None and d.shape[0] > 0:
                    smp = d[:obs_sample]
                    smp = [x.decode() if isinstance(x, bytes) else str(x) for x in smp]
                    print(f"    {c}: {str(d.dtype)} e.g. {smp}")

        # var
        var = f["var"]
        vcols = list(var.keys())
        vid = var.get("_index")
        n_genes = vid.shape[0] if vid is not None else None
        gene_sample = [x.decode() if isinstance(x, bytes) else str(x) for x in vid[:10]] if vid is not None else []
        print(f"  var: n_genes={n_genes} columns={vcols} index e.g. {gene_sample}")
        for c in vcols:
            g = var[c]
            d = g if isinstance(g, h5py.Dataset) else g.get("data")
            if d is not None and d.shape[0] > 0:
                smp = [x.decode() if isinstance(x, bytes) else str(x) for x in d[:5]]
                print(f"    {c}: e.g. {smp}")

        # quick classification of columns
        ident = [c for c in cols if any(h in c.lower() for h in IDENTITY_HINTS)]
        demo = [c for c in cols if any(h in c.lower() for h in DEMO_HINTS)]
        cts = [c for c in cols if any(h in c.lower() for h in CT_HINTS)]
        print(f"  ID-candidate cols: {ident}")
        print(f"  demographic-candidate cols: {demo}")
        print(f"  cell-type-candidate cols: {cts}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    audit(args[0])
