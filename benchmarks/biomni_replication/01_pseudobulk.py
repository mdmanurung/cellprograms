#!/usr/bin/env python
"""Donor x cell-type pseudobulk replication of Biomni's cohorts (COMBAT recipe).

Usage: 01_pseudobulk.py --cohort {combat,stephenson,hlca,onek1k} [--source {X,raw}] [--obs-only]

Dataset is resolved from the pinned CELLxGENE Census datasets table, the source
h5ad is fetched with cellxgene_census.download_source_h5ad (Census obs has no
author labels such as COMBAT major_subset / Stephenson author_cell_type, so the
source h5ad is the only place for them), the matrix is streamed in row chunks and
summed per (donor, cell type) with a sparse indicator matmul.
--source X (default, output data/<cohort>/): sums the log-normalized `X` (log1p CP10k) AS-IS.
  This is what reproduces Biomni's trace: its pb_counts were sums of log-normalized values, and its
  prep_matrix (verbatim below) then applied CPM + log2 to those sums, after truncating the sums to integers
  (np.floor). Verified: reproduces the trace cell-4 marker table exactly (max abs diff 0.00 over 120 CT x marker cells).
--source raw (output data/<cohort>_raw/): sums raw counts (raw.X, or X when the h5ad has no raw, as OneK1K).
--obs-only: obs-level QC (donor counts per CT vs reference) only, no counts.
"""
import argparse, json, os
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = "2025-11-08"
MIN_CELLS, N_HVG, MIN_FRAC = 50, 3000, 0.10
MISSING = {"nan", "NaN", "NA", "None", "", "unknown"}  # string-typed NaN labels in h5ad obs

COMBAT_CTS = ["CD4", "CD8", "B", "NK", "cMono", "ncMono", "DC", "PB", "GDT", "DP"]
COMBAT_GROUPS = {"T": ["CD4", "CD8", "GDT", "DP"], "B_PB": ["B", "PB"], "NK": ["NK"],
                 "Mono": ["cMono", "ncMono"], "DC": ["DC"]}
HLCA_CTS = {  # Census cell_type -> short name used in the reference report
    "alveolar macrophage": "AM", "elicited macrophage": "EMac", "classical monocyte": "cMono",
    "CD1c-positive myeloid dendritic cell": "DC", "CD8-positive, alpha-beta T cell": "CD8",
    "multiciliated columnar cell of tracheobronchial tree": "MCC",
    "pulmonary alveolar type 2 cell": "AT2", "respiratory basal cell": "RBC",
    "vein endothelial cell": "VEC"}

ONEK1K_CTS = {  # Census cell_type -> short name used in the reference report
    "central memory CD4-positive, alpha-beta T cell": "CD4_CM", "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4_Naive",
    "effector memory CD4-positive, alpha-beta T cell": "CD4_EM", "effector memory CD8-positive, alpha-beta T cell": "CD8_EM",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8_Naive", "naive B cell": "B_Naive", "memory B cell": "B_Mem",
    "natural killer cell": "NK", "CD14-positive monocyte": "CD14_Mono"}

# label: obs column; min_donors: keep CTs with >= this many donors at >=50 cells (None = fixed list)
CFG = {
    "combat": dict(
        select=lambda d: d[d.collection_id.eq("8f126edf-5405-4731-8374-b5ce11f53e82")
                           & d.dataset_id.str.startswith("ebc2e1ff")],
        label="major_subset", cts=COMBAT_CTS, groups=COMBAT_GROUPS, agg="first",
        meta=["donor_id", "Source", "disease", "Age", "BMI", "Hospitalstay", "Death28", "Outcome",
              "Institute", "SARSCoV2PCR", "Symptomatic", "TimeSinceOnset", "Smoking",
              "Requiredvasoactive", "sex", "development_stage", "scRNASeq_sample_ID", "GEX_region"]),
    "stephenson": dict(
        select=lambda d: d[d.collection_id.eq("ddfad306-714d-4cc0-9985-d9072820c530")],
        label="author_cell_type", min_donors=30, agg="first",
        meta=["donor_id", "Status", "Worst_Clinical_Status", "Outcome", "Smoker", "sex", "Days_from_onset",
              "development_stage", "Site", "Collection_Day", "Resample", "disease", "sample_id"]),
    "hlca": dict(
        select=lambda d: d[d.collection_id.eq("6f6d381a-7701-4781-935c-db10d30de293")
                           & d.dataset_title.str.contains(r"\(core\)")],
        label="cell_type", cts=list(HLCA_CTS), rename=HLCA_CTS, agg="mode",
        primary=False,  # HLCA core h5ad has is_primary_data==False for all cells (primary copies live in the "full" dataset)
        meta=["donor_id", "tissue_level_2", "smoking_status", "study", "dataset", "sequencing_platform",
              "age_or_mean_of_age_range", "sex", "BMI", "disease", "subject_type"]),
    "onek1k": dict(
        select=lambda d: d[d.collection_id.eq("dde06e0f-ab3b-46be-96a2-a8082383c4a1")],
        label="cell_type", cts=list(ONEK1K_CTS), rename=ONEK1K_CTS, agg="first", cap=300, seed=42,
        meta=["donor_id", "sex", "age", "development_stage", "pool_number", "disease"]),
}


def load_obs(cohort):
    import cellxgene_census
    cfg = CFG[cohort]
    with cellxgene_census.open_soma(census_version=CENSUS) as c:
        ds = c["census_info"]["datasets"].read().concat().to_pandas()
    row = cfg["select"](ds)
    assert len(row) == 1, row[["dataset_id", "dataset_title"]]
    dsid = row.dataset_id.iloc[0]
    path = f"{HERE}/data/_h5ad/{dsid}.h5ad"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        cellxgene_census.download_source_h5ad(dsid, to_path=path + ".part", census_version=CENSUS)
        os.rename(path + ".part", path)
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs.copy()
    print(f"{cohort}: dataset {dsid} {a.shape}", flush=True)
    return a, obs, dsid


def primary(cohort, obs):
    return obs[obs.is_primary_data.astype(bool)] if CFG[cohort].get("primary", True) else obs


def select_samples(cohort, obs, log):
    """Return obs restricted to kept cells with columns donor, ct; plus sample table."""
    cfg = CFG[cohort]
    obs = primary(cohort, obs)
    if "disease" in obs and cohort == "onek1k":
        log["onek1k_disease_counts"] = obs.drop_duplicates("donor_id").disease.astype(str).value_counts().to_dict()
        obs = obs[obs.disease.astype(str) == "normal"]
    if cfg.get("cap"):  # inferred rule, see README
        donors = np.array(sorted(obs.donor_id.astype(str).unique()))
        if len(donors) > cfg["cap"]:
            keep = np.random.default_rng(cfg["seed"]).permutation(donors)[:cfg["cap"]]
            obs = obs[obs.donor_id.astype(str).isin(set(keep))]
        log["cap_rule"] = f"seeded default_rng(seed={cfg['seed']}).permutation(sorted donors)[:{cfg['cap']}] over healthy donors [INFERRED]"
    lab = obs[cfg["label"]].astype(str)
    obs = obs.assign(donor=obs.donor_id.astype(str), ct=lab.where(~lab.isin(MISSING)))
    n = obs.dropna(subset=["ct"]).groupby(["donor", "ct"], observed=True).size()
    n = n[n >= MIN_CELLS]
    if "cts" in cfg:
        keep_ct = cfg["cts"]
    else:
        nd = n.groupby(level="ct").size()
        keep_ct = sorted(nd[nd >= cfg["min_donors"]].index) if cfg["min_donors"] else sorted(nd.index)
    n = n[n.index.get_level_values("ct").isin(keep_ct)]
    return obs, n, keep_ct


def group_sum(a, obs, n, source="X", chunk=100_000):
    """genes x samples matrix via streamed sparse group-sum of X (as-is) or raw.X (counts)."""
    samples = [f"{d}||{c}" for d, c in n.index]
    sid = {k: i for i, k in enumerate(n.index)}
    code = np.full(a.n_obs, -1, dtype=np.int64)
    pos = pd.Series(np.arange(a.n_obs), index=a.obs_names)
    key = list(zip(obs.donor, obs.ct))
    idx = pos[obs.index].to_numpy()
    code[idx] = [sid.get(k, -1) for k in key]
    use_raw = source == "raw" and a.raw is not None
    X = a.raw.X if use_raw else a.X
    var = a.raw.var if use_raw else a.var
    acc = np.zeros((len(samples), var.shape[0]))
    for s in range(0, a.n_obs, chunk):
        cc = code[s:s + chunk]
        m = cc >= 0
        if not m.any():
            continue
        x = X[s:s + chunk]
        x = sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr()
        if s == 0 and source == "raw":
            assert np.allclose(x.data[:10000], np.round(x.data[:10000])), "raw.X not integer counts"
        ind = sp.csr_matrix((np.ones(m.sum()), (cc[m], np.flatnonzero(m))), shape=(len(samples), len(cc)))
        acc += (ind @ x).toarray()
        print(f"  streamed {min(s + chunk, a.n_obs)}/{a.n_obs}", flush=True)
    if source == "X":  # Biomni's pb_counts were these sums cast to integer (truncated); floor reproduces the trace marker table exactly (max diff 0.00)
        acc = np.floor(acc)
    return pd.DataFrame(acc.T, index=var.index, columns=samples), var


def prep_matrix(c, n_hvg=N_HVG, min_frac=MIN_FRAC):  # verbatim COMBAT recipe; c genes x samples
    c = c.loc[(c >= 1).mean(axis=1) >= min_frac]
    lc = np.log2(c.div(c.sum(axis=0), axis=1) * 1e6 + 1).T
    v = lc.var(axis=0).sort_values(ascending=False)
    return lc[v.index[:min(n_hvg, len(v))]]


def donor_meta(cohort, obs_all):
    cfg = CFG[cohort]
    cols = [c for c in cfg["meta"] if c in obs_all.columns]
    o = obs_all[cols].copy()
    for c in o.columns:  # categorical -> object (string "nan" -> NaN)
        o[c] = o[c].astype(object).replace({"nan": np.nan})
    if cfg["agg"] == "first":
        m = o.drop_duplicates("donor_id").set_index("donor_id")
    else:  # mode of non-null values per donor
        m = o.groupby("donor_id", observed=True).agg(
            lambda s: s.dropna().mode().iloc[0] if s.notna().any() else np.nan)
    if cohort == "stephenson":  # Age INFERRED from development_stage "<N>-year-old stage"; NaN for decade-only stages
        m["Age"] = pd.to_numeric(m["development_stage"].astype(str).str.extract(r"^(\d+)-year-old")[0])
    if cohort == "hlca":
        m["Age"] = pd.to_numeric(m["age_or_mean_of_age_range"], errors="coerce")
    return m


def composition(obs, label_col):
    lab = obs[label_col].astype(str).where(lambda s: ~s.isin(MISSING), "nan")  # unlabelled kept as "nan" column
    ct = pd.crosstab(obs.donor_id.astype(str), lab)
    return ct.div(ct.sum(axis=1), axis=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True, choices=list(CFG))
    p.add_argument("--source", choices=["X", "raw"], default="X")
    p.add_argument("--obs-only", action="store_true")
    args = p.parse_args()
    cohort, cfg = args.cohort, CFG[args.cohort]
    out = f"{HERE}/data/{cohort}" + ("_raw" if args.source == "raw" else "")
    os.makedirs(f"{out}/mat", exist_ok=True)

    a, obs, dsid = load_obs(cohort)
    log = {"cohort": cohort, "dataset_id": dsid, "census_version": CENSUS, "source": args.source, "n_cells_total": int(a.n_obs)}
    obs_sel, n, keep_ct = select_samples(cohort, obs, log)
    rename = cfg.get("rename", {})
    short = lambda c: rename.get(c, c)
    donors_per_ct = {short(c): int(n.xs(c, level="ct").size) for c in keep_ct}

    ref = pd.read_csv(f"{HERE}/ref/cellprograms_benchmark/{cohort}/tables/program_summary.csv")
    ref = ref[ref.variant == "none"].set_index("cell_type").n_obs
    qc = {}
    for c in sorted(set(donors_per_ct) | set(ref.index)):
        got, exp = donors_per_ct.get(c), int(ref[c]) if c in ref.index else None
        qc[c] = {"n_donors": got, "ref_n_obs": exp, "status": "OK" if got == exp else "MISMATCH"}
        print(f"{c:20s} got={got} ref={exp} {qc[c]['status']}", flush=True)
    log.update(qc=qc, all_ok=all(v["status"] == "OK" for v in qc.values()),
               n_primary_cells=int(len(obs_sel)), n_donors=int(obs_sel.donor.nunique()))
    json.dump(log, open(f"{out}/qc.json", "w"), indent=1)
    if args.obs_only:
        return

    obs_sel = obs_sel[obs_sel.ct.notna()]
    counts, var = group_sum(a, obs_sel, n, args.source)
    ncell = n.rename(lambda x: x).reset_index()
    ncell["sample"] = ncell.donor + "||" + ncell.ct
    coldata = ncell.rename(columns={0: "n_cells", "ct": "cell_type", "donor": "donor_id"})
    coldata["cell_type"] = coldata.cell_type.map(short)
    counts.to_parquet(f"{out}/pb_counts.parquet")
    coldata.to_csv(f"{out}/pb_coldata.csv", index=False)
    vv = var.copy()
    if "feature_name" not in vv:
        vv["feature_name"] = a.var.feature_name.reindex(vv.index)
    vv.rename_axis("feature_id").reset_index().to_csv(f"{out}/pb_var.csv", index=False)

    donor_of = coldata.set_index("sample").donor_id
    ct_of = coldata.set_index("sample").cell_type
    cols_of = lambda cts: [s for s in counts.columns if ct_of[s] in cts]
    for ct in sorted(set(ct_of)):
        cols = cols_of([ct])
        m = prep_matrix(counts[cols]); m.index = donor_of[cols].values
        m.to_csv(f"{out}/mat/ct_{ct}.csv")
        print(f"{ct}: {m.shape}", flush=True)
    glob = counts.T.groupby(donor_of[counts.columns].values).sum().T
    prep_matrix(glob).to_csv(f"{out}/mat/global.csv")
    for gn, members in (cfg.get("groups") or {}).items():
        cols = cols_of(members)
        g = counts[cols].T.groupby(donor_of[cols].values).sum().T
        prep_matrix(g).to_csv(f"{out}/mat/group_{gn}.csv")

    obs_p = primary(cohort, obs)
    if cfg.get("cap"):  # restrict to the capped (healthy, seeded-sampled) donors, not all primary-cell donors
        obs_p = obs_p[obs_p.donor_id.astype(str).isin(set(obs_sel.donor))]
    donor_meta(cohort, obs_p).to_csv(f"{out}/donor_meta.csv")
    composition(obs_p, cfg["label"]).to_csv(f"{out}/composition.csv")
    print("done", flush=True)


if __name__ == "__main__":
    main()
