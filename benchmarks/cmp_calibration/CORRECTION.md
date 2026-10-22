# Correction to the historical confirmation analysis

`CONFIRM_ANALYSIS.md`, CSV summaries and original preregistrations are historical,
not current calibration certificates. The frozen base retains the as-run code.
The old analyzer is retired rather than allowed to regenerate invalid verdicts.

- `any_sig` includes true-effect calls in mixed-null scenarios and is not an
  any-false event. The B1 row's value 1.000 does not estimate false-family error.
- Pooled null-node FPR is not FWER; the “STRONG FWER CONTROL” verdict is withdrawn.
- Node max-test adjustment does not calibrate links or merged CMP memberships.
- Historical family events are unknown unless actual replicate-level true-null
  counts or node IDs can be recovered and their generator/null definition checked.
  Never reconstruct those events from rounded aggregate rates.

`analyze_corrected.R` accepts recovered node counts and a complete explicit
requested grid. It reports M06 and M07 separately, includes missing/failed/
uncallable jobs, and supplies a one-sided exact binomial upper bound for
replicate events, not a binomial node-level CI for dependent nodes. Its historical
known-score scope is conditional on trusted truth/count records; it never
certifies selected-model, edge or set calibration. New campaign adapters should
retain exact true-null node IDs and callability via the reproducibility workflow.
