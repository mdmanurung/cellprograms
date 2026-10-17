"""Workstream F (PILOT): cross-arm program comparison.

Arms (9): {ebmf, cnmf, genenmf} x {pseudobulk, metacell, singlecell}.
Program representation: genes x K loading matrix per arm.
  ebmf    : F_ebmf_{arm}_loadings.csv (genes x K)
  cnmf    : F_cnmf_{arm}_spectra.csv (K x genes) -> transposed
  genenmf : MP-level loadings = sign-aligned mean of member patterns' w
            (members read from mp_composition.csv; fallback = per-K solutions)
Metrics:
  - pairwise Hungarian-matched |cor| of gene loadings (36 arm pairs) -> heatmap
  - top-50 gene-set Jaccard on matched pairs
  - marker-set rank enrichment per program
  - donor-score concordance (aligned Pearson) + age association (Spearman)
"""
import itertools
import json
import os

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.optimize import linear_sum_assignment

ROOT = "/mnt/shared-workspace/shared/audit/F_granularity"
ARMS = ["pseudobulk", "metacell", "singlecell"]
METHODS = ["ebmf", "cnmf", "genenmf", "cebmf"]
TOPN = 50

MARKERS = {
    "Cytotoxicity": ["GZMB", "GZMK", "PRF1", "GZMA", "GZMH", "NKG7", "GNLY", "FGFBP2", "KLRD1", "KLRK1"],
    "IFN_response": ["STAT1", "IRF7", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "IFI6", "IFI44L"],
    "Cell_cycle": ["MKI67", "TOP2A", "STMN1", "TYMS", "PTTG1", "AURKB"],
    "Naive_memory": ["TCF7", "IL7R", "LEF1", "CCR7", "SELL"],
    "Terminal_effector": ["CX3CR1", "ZEB2", "ID2", "S1PR5"],
    "Heat_shock": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPB1", "HSPH1"],
    "Antigen_presentation": ["HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1"],
    "NFkB_instant": ["NFKBIA", "TNFAIP3", "FOS", "JUN", "JUNB"],
}


def hungarian_match(W1, W2):
    """Match programs by |cor| of gene loadings. Returns list of (i, j, cor_signed, cor_abs)."""
    C = np.corrcoef(W1.T, W2.T)[:W1.shape[1], W1.shape[1]:]
    C = np.nan_to_num(C, nan=0.0)
    ri, ci = linear_sum_assignment(-np.abs(C))
    return [(int(i), int(j), float(C[i, j]), float(abs(C[i, j]))) for i, j in zip(ri, ci)]


def top_jaccard(W1, j1, W2, j2, top=TOPN):
    g1 = set(W1.index[np.argsort(-np.abs(W1.iloc[:, j1].values))[:top]])
    g2 = set(W2.index[np.argsort(-np.abs(W2.iloc[:, j2].values))[:top]])
    return len(g1 & g2) / len(g1 | g2)


def marker_enrichment(W, markers):
    """Mean rank percentile of marker genes per program (1 = top of |loading|)."""
    ranks = W.abs().rank(pct=True)  # genes x K
    out = {}
    for name, genes in markers.items():
        gs = [g for g in genes if g in ranks.index]
        out[name] = float(ranks.loc[gs].mean().mean()) if gs else np.nan
    return out


def load_arm(method, arm):
    """Returns (W genes x K, donor_scores DataFrame donors x K, K)."""
    if method in ("ebmf", "cebmf"):
        W = pd.read_csv(f"{ROOT}/F_{method}_{arm}_loadings.csv", index_col=0)
        D = pd.read_csv(f"{ROOT}/F_{method}_{arm}_donor_scores.csv", index_col=0)
        D.index = D.index.astype(str)
        return W, D, W.shape[1]
    if method == "cnmf":
        spc = pd.read_csv(f"{ROOT}/F_cnmf_{arm}_spectra.csv", index_col=0)
        W = spc.T
        usage = pd.read_csv(f"{ROOT}/F_cnmf_{arm}_usage.csv", index_col=0)
        meta = pd.read_csv(f"{ROOT}/gran/{arm}/obs_meta.tsv", sep="\t", index_col=0)
        meta.index = meta["id"].astype(str).values
        donor = meta.loc[usage.index.astype(str), "donor"].astype(str)
        D = usage.groupby(donor.values).mean()
        return W, D, W.shape[1]
    # genenmf: MP-level W/H from memberships
    # allW: genes x all_patterns; allH: all_patterns x obs (rbind-stacked)
    allW = pd.read_csv(f"{ROOT}/F_genenmf_{arm}_allW.csv", index_col=0)
    allH = pd.read_csv(f"{ROOT}/F_genenmf_{arm}_allH.csv", index_col=0)
    comp = pd.read_csv(f"{ROOT}/F_genenmf_{arm}_mp_composition.csv")
    comp.columns = [c.strip() for c in comp.columns]
    patcol = next(c for c in comp.columns if comp[c].astype(str).str.contains("pattern").any())
    mpcol = next(c for c in comp.columns if c.lower().startswith("mp") or "metaprogram" in c.lower())
    meta = pd.read_csv(f"{ROOT}/gran/{arm}/obs_meta.tsv", sep="\t", index_col=0)
    meta.index = meta["id"].astype(str).values
    donor = meta["donor"].astype(str)
    Ws, Ds, kept = [], [], []
    for mp_name, grp in comp.groupby(comp[mpcol].astype(str)):
        pats = [p for p in grp[patcol].astype(str) if p in allW.columns]
        if not pats:
            continue
        Wm = allW[pats].copy()          # genes x members
        Hm = allH.loc[pats].copy()      # members x obs
        signs = [1.0]
        for c in pats[1:]:
            r = np.corrcoef(Wm[pats[0]].values, Wm[c].values)[0, 1]
            s = -1.0 if r < 0 else 1.0
            signs.append(s)
            Wm[c] = Wm[c] * s
            Hm.loc[c] = Hm.loc[c] * s
        Ws.append(Wm.mean(axis=1))
        Ds.append(Hm.mean(axis=0))      # per-obs MP usage
        kept.append(mp_name)
    W = pd.DataFrame(Ws).T
    W.columns = kept
    D = pd.DataFrame({mp: pd.Series(np.asarray(d), index=[str(x) for x in d.index])
                      for mp, d in zip(kept, Ds)})
    D = D.groupby([donor.get(i, i) for i in D.index]).mean()
    return W, D, W.shape[1]


def main():
    recipe = json.load(open(f"{ROOT}/gran/recipe.json"))
    hvg = recipe["hvg"]
    arms_data = {}
    for m in METHODS:
        for a in ARMS:
            try:
                W, D, K = load_arm(m, a)
                genes_common = W.index.intersection(hvg)
                arms_data[(m, a)] = (W.loc[genes_common], D, K)
                print(f"loaded {m}/{a}: K={K} genes={W.shape[0]}")
            except Exception as e:
                print(f"FAILED to load {m}/{a}: {e}")

    names = list(arms_data)
    # pairwise concordance
    rows = []
    for (n1, n2) in itertools.combinations(names, 2):
        W1, D1, _ = arms_data[n1]
        W2, D2, _ = arms_data[n2]
        g = W1.index.intersection(W2.index)
        matches = hungarian_match(W1.loc[g], W2.loc[g])
        cors = [c for _, _, _, ca in matches for c in [ca]]
        jacs = [top_jaccard(W1.loc[g], i, W2.loc[g], j) for i, j, _, _ in matches]
        rows.append({"arm1": f"{n1[0]}_{n1[1]}", "arm2": f"{n2[0]}_{n2[1]}",
                     "K1": W1.shape[1], "K2": W2.shape[1],
                     "mean_matched_abs_cor": float(np.mean(cors)),
                     "min_matched_abs_cor": float(np.min(cors)),
                     "mean_top50_jaccard": float(np.mean(jacs))})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(f"{ROOT}/F_comparison_pairs.csv", index=False)

    # concordance heatmap (9x9)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
    matplotlib.rcParams["svg.fonttype"] = "none"
    M = np.eye(len(names))
    lookup = {(f"{n1[0]}_{n1[1]}", f"{n2[0]}_{n2[1]}"): v for n1, n2, v in
              [(tuple(r[["arm1", "arm2"]]), r["mean_matched_abs_cor"]) for _, r in pairs.iterrows()]}
    for i, n1 in enumerate(names):
        for j, n2 in enumerate(names):
            if i == j:
                M[i, j] = 1.0
            else:
                key = (f"{n1[0]}_{n1[1]}", f"{n2[0]}_{n2[1]}")
                key_r = (key[1], key[0])
                M[i, j] = lookup.get(key, lookup.get(key_r, np.nan))
    labels = [f"{m}_{a}" for m, a in names]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(M, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(names)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(names)), labels)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                    color="white" if M[i, j] < 0.6 else "black", fontsize=8)
    fig.colorbar(im, label="mean matched |cor| (gene loadings)")
    fig.tight_layout()
    fig.savefig(f"{ROOT}/F_concordance_heatmap.png", dpi=200)
    plt.close(fig)

    # marker enrichment per arm
    mrows = []
    for (m, a), (W, D, K) in arms_data.items():
        enr = marker_enrichment(W, MARKERS)
        for k in range(K):
            row = {"method": m, "arm": a, "program": f"{m}_{a}_P{k + 1}"}
            row.update({name: float(W.abs().rank(pct=True).loc[
                [g for g in genes if g in W.index], k + 1].mean()) if any(g in W.index for g in genes) else np.nan
                for name, genes in MARKERS.items()})
            mrows.append(row)
    pd.DataFrame(mrows).to_csv(f"{ROOT}/F_marker_enrichment.csv", index=False)

    # donor-score concordance + age
    age = pd.read_csv(f"{ROOT}/gran/donor_age.tsv", sep="\t", index_col=0).iloc[:, 0]
    arows = []
    for (m, a), (W, D, K) in arms_data.items():
        common = D.index.intersection(age.index)
        if len(common) < 10:
            continue
        from scipy.stats import spearmanr
        for k in range(D.shape[1]):
            rho, p = spearmanr(D.loc[common, k], age.loc[common])
            arows.append({"method": m, "arm": a, "program": k + 1,
                          "spearman_age": float(rho), "p": float(p), "n_donors": len(common)})
    pd.DataFrame(arows).to_csv(f"{ROOT}/F_age_assoc.csv", index=False)
    print("comparison done")


if __name__ == "__main__":
    main()
