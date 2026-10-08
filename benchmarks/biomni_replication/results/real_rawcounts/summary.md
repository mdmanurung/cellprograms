# Real COMBAT: PA-scores vs multi-view PMD (repo K=10 fit)
Script `06_real_pa_vs_pmd.R` (job 25728237). 10 CTs x 10 programs, 124 donors (49-121 per CT).
PA: n_perm=200. PMD: 5 MCPs, n_perm=100, seed 43. "Resid" = each score column residualized on sex + Institute
within the CT's own donors (no NA covariates exist, so the available-case fallback never fired). Source NOT residualized.
Pair rules copied from 05: PA shared iff n_shared_05>=1; PMD shared iff an MCP has p<0.05 and pair_cor>0.5.
Files: pa_pairs_{raw,resid}.csv, pmd_mcps_{raw,resid}.csv, pmd_mcp_assoc_raw_vs_resid.csv, pmd_mcp_signature_hits_*.csv, pmd_pairs_*.csv.

## Pair calls (45 CT pairs)
| arm | PA shared | PMD shared |
|---|---|---|
| raw | 39/45 | 45/45 |
| resid | 38/45 | 45/45 |
- PA: 6 pairs are uncallable in BOTH arms: DC-DP, DC-GDT, DC-PB, DP-GDT, DP-PB, GDT-PB. They have only 18-39 common
  donors, so the random-subspace reference (sqrt(ka)+sqrt(kb))/sqrt(n) is >1 and no cosine can pass it. Among the 39
  callable pairs, PA calls every one shared raw and 38 after residualization (only ncMono-PB lost; its reference is 0.95).
- Mean shared directions per pair 4.2 -> 3.7; mean smallest principal angle 8.7 -> 11.2 deg.
- PMD called all 45 pairs in both arms. All 5 MCPs hit the permutation floor (p=0.0099 = 1/101); mean pairwise r among
  active CTs 0.73-0.93 raw, 0.52-0.91 resid.
- Reading: on real data neither pair call discriminates. PA is limited by donor overlap, PMD saturates. Residualizing
  sex + Institute barely changes pair-level sharing; it changes what the PMD axes are.

## PMD MCPs (state = mean of per-CT MCP scores); p-values uncorrected
Raw:
| MCP | top genes (top-3 CTs) | sex p | Inst p | Source p | Outcome p | Age rho |
|---|---|---|---|---|---|---|
| 01 | Ig V genes (IGLV2-14, IGKV3D-20; DC, CD8, CD4) | 0.10 | 0.14 | 0.76 | 0.83 | 0.01 |
| 02 | Y genes + XIST (DC, CD8, cMono) | 2e-18 | 4e-4 | 2e-6 | 3e-5 | 0.13 |
| 03 | IFI44L IFIT1/3 RSAD2 MX1 (CD8, NK, ncMono) | 0.50 | 0.37 | 2e-3 | 9e-4 | -0.16 |
| 04 | JUN FOS NR4A2/3 TNFAIP3 RGS1 CISH; ncMono GSTM1 | 0.13 | 9e-3 | 2e-16 | 2e-8 | -0.21 |
| 05 | Ig V genes (IGHV1-69D, IGLV3-10) + GSTM1 in DC | 0.89 | 0.08 | 0.30 | 0.46 | -0.03 |
Resid:
| MCP | top genes | sex p | Inst p | Source p | Outcome p | Age rho |
|---|---|---|---|---|---|---|
| 01 | IFN genes (same as raw 03) | 0.95 | 0.68 | 3e-3 | 3e-3 | -0.14 |
| 02 | Ig V genes (same as raw 01) | 0.73 | 0.40 | 0.92 | 0.89 | -0.02 |
| 03 | Y genes + JUN, RGS1, TENT5C, CISH | 0.89 | 0.77 | 2e-15 | 1e-9 | -0.27 |
| 04 | Ig V genes + GSTM1 in DC (same as raw 05) | 0.97 | 0.76 | 0.81 | 0.74 | -0.01 |
| 05 | cell cycle (RRM2 BIRC5 CCNB2 PLK1; CD4); Y genes in CD8/cMono | 0.89 | 0.75 | 5e-12 | 2e-10 | 0.26 |

## Biomni archetypes in our K=10 fit (matched by top genes and associations)
- MCP01 Ig: present raw (our MCP01) and survives (resid MCP02). It had no sex, Institute or Source association, so
  residualizing could not remove it. Ig V genes in 9 CTs point to ambient/plasma-cell RNA or germline variation, not biology of interest.
- MCP05 GSTM1: present raw (our MCP05; GSTM1 in DC, ncMono GSTM1 also in MCP04) and survives (resid MCP04). Same reason;
  its Ig-V and GSTM1 content looks genotype-like (GSTM1 deletion, IGHV copy number), which is a guess we did not test.
- MCP02 sex: present raw (our MCP02, sex p 2e-18). After residualization no MCP associates with sex (all p>=0.73), so the
  sex-association vanishes. Y genes still appear in the traced loadings of resid MCP03/05, because those programs still
  carry within-sex variation. Gene trace shows loadings, not state-sex association.
- MCP04 immediate-early/Institute: present raw (our MCP04). Its Institute association drops (p 9e-3 -> 0.77), but a
  JUN/RGS1/TENT5C/CISH axis stays, bundled with Y-gene programs (resid MCP03), and stays strongly Source/Outcome
  associated (Source p 2e-15). So Institute was only a minor part of it. It does not "vanish"; most of it tracks disease Source.
- Interferon: survives (raw MCP03 -> resid MCP01), same genes in all 10 CTs, Source p 2e-3 and 3e-3. It is the only
  MCP whose Source/Outcome association is modest and stable. A new proliferation MCP (resid MCP05) shows a stronger Source link.

## Caveats and nulls
- Null: PA and PMD both call (almost) every pair shared; all 5 PMD p-values sit at the 1/101 floor. No evidence here that either separates real sharing
  from broad donor-level coordination.
- The sex+Institute removal is per CT on that CT's donors. Sex and Institute are confounded with Source (MCP02 Source
  p 2e-6), so some Source signal is removed along with them.
- fork trace_state() has no 'pmd' branch. We called it through an ICA-shaped integ holding PMD program weights divided
  by each program's score SD (PMD runs on standardized scores). The fork code was not edited.
