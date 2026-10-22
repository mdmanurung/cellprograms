
## 2026-10-22 — MR-2 cNMF arm input + API (see PREREG_mr2-cnmf-v1.md Addendum A)
- cNMF consumes planted counts (round+clip ≥0), not the centered/residualized K30 matrix;
  cNMF requires non-negative integer counts. Each backend tested on its natural input.
- cnmf_fit.py rewritten for cnmf 1.7.1 (silhouette stability, load_results API);
  CNMF_SCRIPT HERE-detection NA guard; writeMM double-extension fix.
- semi_base contract completed: ncell / site / msd fields added (were absent; caused
  "attempt to set an attribute on NULL" and downstream gen_semi failures).
