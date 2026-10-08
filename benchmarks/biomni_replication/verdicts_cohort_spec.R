## verdicts_cohort_spec.R -- explicit per-cohort covariate roles for summarize_eval.R.
## Derived from the covariate names in results/eval/<cohort>_*_corrected.csv and
## ref/cellprograms_benchmark/<cohort>/report_benchmark_<cohort>.md.
##   bio_expr : biology, expression-driven (cp expected to be able to win)
##   bio_comp : biology, compositional by design (severity / disease / smoking / mortality axes; the
##              reports show composition/clr win these, so a cp loss here is not a program failure)
##   bio_other: biology, neither clearly (tallied under "other")
##   tech     : technical/batch (Panel B, lower is better; never in the biology tally)
##   excl     : named vector covariate -> reason (degenerate or duplicate; never tallied)
## A covariate present in the CSVs but absent from all lists is reported as "unclassified" with a warning.
COHORT_SPEC <- list(
  combat = list(
    bio_expr = c("sex", "Age_num", "TimeSinceOnset"),           # TimeSinceOnset/Hospitalstay: judgement call
    bio_comp = c("Source", "Outcome", "Death28"),                # report: Death28 near-purely compositional
    bio_other = c("Hospitalstay"),
    tech = c("Institute", "GEX_region"), excl = c()),
  hlca = list(
    bio_expr = c("sex", "Age", "BMI"),
    bio_comp = c("smoking_status", "subject_type"),             # report: smoking is cell-compositional
    bio_other = c("tissue_level_2"),                             # study-confounded in HLCA
    tech = c("study", "sequencing_platform"),
    excl = c(age_or_mean_of_age_range = "duplicate of Age (identical values)")),
  onek1k = list(
    bio_expr = c("sex", "age"),
    bio_comp = c(), bio_other = c(),
    tech = c("pool_number"), excl = c()),                        # numeric technical covariate
  stephenson = list(
    bio_expr = c("sex", "Age"),
    bio_comp = c("Status", "disease", "Worst_Clinical_Status", "Outcome", "Smoker"),
    bio_other = c(),
    tech = c("Site"),
    excl = c(Collection_Day = "degenerate (119/120 donors one class; every representation = floor)",
             Resample = "degenerate (119/120 donors one class; every representation = floor)"))
)
