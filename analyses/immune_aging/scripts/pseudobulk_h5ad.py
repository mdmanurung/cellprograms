#!/usr/bin/env python3
"""Pseudobulk a portal *_reannotated.h5ad to donor x cell-type raw-count matrices.

Reads .raw/X (CSR) in contiguous row chunks, aggregates sum of counts per
(donor, cell type) for ALL cells, then drops 'Bad quality' groups from output and
from donor cell counts, applies a min-cells-per-group QC filter, and writes:

  <outdir>/pseudobulk.h5      group per cell type: /<safe_ct>/counts (donors x genes
                              float32), /<safe_ct>/donors (strings); plus /genes
  <outdir>/sample_metadata.csv  donor-level: donor_id, age, sex, ethnicity, CMV, BMI,
                              cohort, n_cells (kept cells only)
  <outdir>/ct_summary.csv       cell type, n_donors, n_cells, min/median cells per donor

Usage: python3 pseudobulk_h5ad.py <in.h5ad> <outdir> <cohort> [--min-cells 10]
"""
import argparse
import csv
import os
import re

import h5py
import numpy as np
import scipy.sparse as sp


def read_categorical(f, name):
    g = f[f"obs/{name}"]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]
        codes = g["codes"][:]
        return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
    d = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in d], dtype=object)


def read_numeric(f, name):
    g = f[f"obs/{name}"]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]
        return np.array([cats[c] for c in g["codes"][:]], dtype=object)
    return g[:]


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")[:80]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("cohort")
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--ct-col", default="annotation")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--chunk", type=int, default=25_000)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    f = h5py.File(args.input, "r")

    donor = read_categorical(f, args.donor_col)
    ct = read_categorical(f, args.ct_col)
    age = read_numeric(f, "age")
    sex = read_categorical(f, "sex")
    eth = read_categorical(f, "ethnicity")
    cmv = read_categorical(f, "CMV")
    bmi = read_categorical(f, "BMI")

    keep_cell = np.array([not c.lower().startswith("bad quality") for c in ct])
    print(f"cells total={len(ct)}, kept (non bad-quality)={int(keep_cell.sum())}")

    # group ids over ALL cells (contiguous chunks); bad-quality groups dropped at write
    pairs = {}
    pair_index = np.empty(len(ct), dtype=np.int64)
    for i, (d, c) in enumerate(zip(donor, ct)):
        key = (d, c)
        idx = pairs.get(key)
        if idx is None:
            idx = len(pairs)
            pairs[key] = idx
        pair_index[i] = idx
    n_groups = len(pairs)
    group_keys = sorted(pairs, key=pairs.get)
    donors_all = list(dict.fromkeys(d for d, _ in group_keys))
    donor_pos = {d: i for i, d in enumerate(donors_all)}
    cts = list(dict.fromkeys(c for _, c in group_keys))
    ct_is_bad = {c: c.lower().startswith("bad quality") for c in cts}
    print(f"groups={n_groups}, donors={len(donors_all)}, cell_types={len(cts)} "
          f"({sum(ct_is_bad.values())} bad-quality types dropped at write)")

    genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/_index"][:]]
    n_genes = len(genes)

    X = f["raw/X"]
    indptr = X["indptr"][:]
    n_cells_h5 = len(indptr) - 1
    assert n_cells_h5 == len(ct), f"row mismatch: {n_cells_h5} vs obs {len(ct)}"

    # per-ct accumulation arrays: donor x genes (skip bad-quality types entirely)
    ct_rows = {c: np.array([donor_pos[d] for d, cc in group_keys if cc == c])
               for c in cts if not ct_is_bad[c]}
    acc = {c: np.zeros((len(ix), n_genes), dtype=np.float32)
           for c, ix in ct_rows.items()}
    group_ct_row = {}
    for gi, (d, c) in enumerate(group_keys):
        group_ct_row[gi] = None if ct_is_bad[c] else (c, int(np.where(ct_rows[c] == donor_pos[d])[0][0]))

    n_cells = len(ct)
    for s in range(0, n_cells, args.chunk):
        e = min(s + args.chunk, n_cells)
        lo, hi = int(indptr[s]), int(indptr[e])
        data = X["data"][lo:hi]
        indices = X["indices"][lo:hi]
        ptr = indptr[s:e + 1] - lo
        xb = sp.csr_matrix((data, indices, ptr), shape=(e - s, n_genes), dtype=np.float64)
        gids = pair_index[s:e]
        uniq, inv = np.unique(gids, return_inverse=True)
        G = sp.csr_matrix((np.ones(len(inv), dtype=np.float64), inv,
                           np.arange(len(inv) + 1)), shape=(len(inv), len(uniq)))
        agg = (G.T @ xb).toarray()
        for j, gi in enumerate(uniq):
            tgt = group_ct_row[int(gi)]
            if tgt is not None:
                c, r = tgt
                acc[c][r] += agg[j]
        print(f"  chunk {e}/{n_cells} done", flush=True)
    f.close()

    # write h5 (float32)
    h5path = os.path.join(args.outdir, "pseudobulk.h5")
    with h5py.File(h5path, "w") as h:
        dt = h5py.string_dtype(encoding="utf-8")
        h.create_dataset("genes", data=np.array(genes, dtype=object), dtype=dt)
        h.create_dataset("donors_all", data=np.array(donors_all, dtype=object), dtype=dt)
        h.create_dataset("cell_types", data=np.array(
            [c for c in cts if not ct_is_bad[c]], dtype=object), dtype=dt)
        for c, ix in ct_rows.items():
            g = h.create_group(safe_name(c))
            g.create_dataset("counts", data=acc[c].astype(np.float32))
            g.create_dataset("donors", data=np.array(
                [donors_all[i] for i in ix], dtype=object), dtype=dt)

    # donor metadata (kept cells only) + per-group cell counts for QC
    cell_counts = np.bincount(pair_index, minlength=n_groups)
    donor_cells = {d: 0 for d in donors_all}
    for (d, c), n in zip(group_keys, cell_counts):
        if not ct_is_bad[c]:
            donor_cells[d] += int(n)

    first_row = {}
    for i, d in enumerate(donor):
        if d not in first_row:
            first_row[d] = i
    meta_rows = []
    for d in donors_all:
        i = first_row[d]
        meta_rows.append({
            "donor_id": d, "age": age[i], "sex": sex[i], "ethnicity": eth[i],
            "CMV": cmv[i], "BMI": bmi[i], "cohort": args.cohort,
            "n_cells": donor_cells[d]})

    with open(os.path.join(args.outdir, "sample_metadata.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(meta_rows[0]))
        w.writeheader()
        w.writerows(meta_rows)

    with open(os.path.join(args.outdir, "ct_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell_type", "n_donors", "n_cells", "min_cells_per_donor",
                    "median_cells_per_donor", "n_donors_pass_qc"])
        for c in cts:
            if ct_is_bad[c]:
                continue
            ns = sorted(int(n) for (d, cc), n in zip(group_keys, cell_counts) if cc == c)
            w.writerow([c, len(ns), sum(ns), ns[0], int(np.median(ns)),
                        sum(1 for n in ns if n >= args.min_cells)])

    print(f"wrote {h5path} and metadata for {len(donors_all)} donors, "
          f"{sum(not ct_is_bad[c] for c in cts)} cell types")


if __name__ == "__main__":
    main()
