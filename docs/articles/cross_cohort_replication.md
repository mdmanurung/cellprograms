# Cross-cohort replication: a preregistered-style case study

This article is an end-to-end, fully evaluated template for asking:
**does a program-level biological effect found in one cohort replicate
in independent cohorts?** It uses a simulated three-cohort fixture so
the truth is known, and it orchestrates existing exports only — no
workflow-specific wrapper. The two transfer pathways are kept distinct:

- **Fixed-basis pathway** — the discovery fit’s learned transform (gene
  panel, centering/scaling, loading effects) is *applied* to the query
  cohort via
  [`project_programs_to_cells()`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md).
  [`freeze_preprocessing()`](https://mdmanurung.github.io/cellprograms/reference/freeze_preprocessing.md)
  alone is **not** this: it freezes decisions (gene panel, design), not
  the learned transform.
- **De-novo pathway** — the query cohort is fit independently (here with
  the frozen feature panel, so the gene universe matches) and programs
  are aligned afterwards with
  [`match_programs()`](https://mdmanurung.github.io/cellprograms/reference/match_programs.md).
  Factor numbers and raw scores are never compared directly across fits.

Technical transportability (coverage, reconstruction) is reported
separately from biological effect replication, and every abstention is
counted with its denominator.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellprograms``)`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``20261019``)`

## Three cohorts, known truth

Two cell types (`T`, `B`), 24 donors per cohort, 400 genes. Planted
structure in T cells:

| Program | Loadings | Activity | Truth |
|----|----|----|----|
| P1 | genes 1-40, identical across cohorts | correlated with donor outcome `y` in every cohort | conserved biological signal |
| P2 | genes 61-100, identical across cohorts | driven by a batch covariate in cohort B only | nuisance-only effect |
| P3 | genes 121-160 | active in cohort C only | cohort-specific program |

B cells carry one conserved program (genes 201-240) with no outcome
link.

\
`n_genes`` ``<-`` ``400``; ``n_d`` ``<-`` ``24`\
`genes`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"G"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_genes``)``)`\
`mk_w`` ``<-`` ``function``(``idx``)`` ``{`` ``w`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``0``, ``n_genes``)``; ``w``[``idx``]`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``idx``)``)``; ``w`` ``}`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`w1`` ``<-`` ``mk_w``(``1``:``40``)``; ``w2`` ``<-`` ``mk_w``(``61``:``100``)``; ``w3`` ``<-`` ``mk_w``(``121``:``160``)``; ``w4`` ``<-`` ``mk_w``(``201``:``240``)`\
\
`mk_cohort`` ``<-`` ``function``(``name``, ``with_p3`` ``=`` ``FALSE``, ``p2_batch`` ``=`` ``FALSE``)`` ``{`\
`  ``donors`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``name``, ``"_d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_d``)``)`\
`  ``z1`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``)``; ``z2`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``)``; ``z4`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``)`\
`  ``batch`` ``<-`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"b1"``, ``"b2"``)``, each ``=`` ``n_d`` ``/`` ``2``)``)`\
`  ``if`` ``(``p2_batch``)`` ``z2`` ``<-`` ``1.2`` ``*`` ``(``batch`` ``==`` ``"b2"``)`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``, sd ``=`` ``0.3``)`\
`  ``y`` ``<-`` ``0.7`` ``*`` ``z1`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``, sd ``=`` ``0.8``)``          ``# outcome linked to P1`\
`  ``YT`` ``<-`` `[`tcrossprod`](https://rdrr.io/r/base/crossprod.html)`(``z1``, ``w1``)`` ``+`` `[`tcrossprod`](https://rdrr.io/r/base/crossprod.html)`(``z2``, ``w2``)`\
`  ``if`` ``(``with_p3``)`` ``YT`` ``<-`` ``YT`` ``+`` `[`tcrossprod`](https://rdrr.io/r/base/crossprod.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d``, sd ``=`` ``1.5``)``, ``w3``)`\
`  ``YT`` ``<-`` ``YT`` ``+`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d`` ``*`` ``n_genes``, sd ``=`` ``0.8``)``, ``n_d``, ``n_genes``)`\
`  ``YB`` ``<-`` `[`tcrossprod`](https://rdrr.io/r/base/crossprod.html)`(``z4``, ``w4``)`` ``+`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_d`` ``*`` ``n_genes``, sd ``=`` ``0.8``)``, ``n_d``, ``n_genes``)`\
`  `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``YT``)`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``YB``)`` ``<-`` ``genes`\
`  `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``YT``)`` ``<-`` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``YB``)`` ``<-`` ``donors`\
`  ``meta`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``observation_id ``=`` ``donors``, cohort ``=`` ``name``, outcome ``=`` ``y``,`\
`                     batch ``=`` ``batch``, stringsAsFactors ``=`` ``FALSE``)`\
`  `[`list`](https://rdrr.io/r/base/list.html)`(``data ``=`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(`[`list`](https://rdrr.io/r/base/list.html)`(``T ``=`` ``YT``, B ``=`` ``YB``)``, sample_metadata ``=`` ``meta``)``,`\
`       meta ``=`` ``meta``)`\
`}`\
`cohA`` ``<-`` ``mk_cohort``(``"A"``)``; ``cohB`` ``<-`` ``mk_cohort``(``"B"``, p2_batch ``=`` ``TRUE``)`\
`cohC`` ``<-`` ``mk_cohort``(``"C"``, with_p3 ``=`` ``TRUE``)`

The predeclared hypothesis, fixed before looking at any query cohort:
**the T-cell program whose activity tracks `outcome` in the discovery
cohort replicates with the same sign as the discovery effect at p \<
0.05.** (A fitted program’s sign is arbitrary; canonicalization fixes it
within a fit, so the rule is same-sign-as-discovery, not a fixed
direction.) Everything else is descriptive.

## Discovery fit (cohort A)

\
`fitA`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``cohA``$``data``, max_factors ``=`` ``4``, seed ``=`` ``1``,`\
`                              feature_scale ``=`` ``"simulated"``)`\
`fitA`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fitA``)`\
[`program_summary`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)`(``fitA``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"cell_type"``, ``"n_programs"``, ``"R2"``)``]`\
`#>   cell_type n_programs         R2`\
`#> 1         T          2 0.26394663`\
`#> 2         B          1 0.07886538`

Freeze the discovery *decisions* (feature panel, centering) for later
de-novo refits, and keep the fitted object itself for fixed-basis
projection:

\
`frozen`` ``<-`` `[`freeze_preprocessing`](https://mdmanurung.github.io/cellprograms/reference/freeze_preprocessing.md)`(``cohA``$``data``, features ``=`` ``"all"``, center ``=`` ``TRUE``)`\
`frozen``$``digest`\
`#> [1] "a92378907d03ffd308f4d967f75bf62bc71c8173e89c91233e00df622462a720"`

## Pathway 1: fixed-basis projection into cohorts B and C

The discovery programs are applied to each query cohort’s donor
pseudobulks. Rows are `type_donor` so both cell types project in one
call; the discovery preprocessing is replayed from the fit, not refit.

\
`stack`` ``<-`` ``function``(``coh``)`` ``{`\
`  `[`rbind`](https://rdrr.io/r/base/cbind.html)`(``coh``$``data``$``pseudobulk``$``T``, ``coh``$``data``$``pseudobulk``$``B``)`` ``|>`\
`    ``` `rownames<-` ```(`[`c`](https://rdrr.io/r/base/c.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(``"T_"``, ``coh``$``meta``$``observation_id``)``,`\
`                   `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"B_"``, ``coh``$``meta``$``observation_id``)``)``)`\
`}`\
`project_cohort`` ``<-`` ``function``(``fit``, ``coh``)`` ``{`\
`  ``expr`` ``<-`` `[`stack`](https://rdrr.io/r/utils/stack.html)`(``coh``)`\
`  ``cts`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`sub`](https://rdrr.io/r/base/grep.html)`(``"_.*"``, ``""``, `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``expr``)``)``, `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``expr``)``)`\
`  `[`project_programs_to_cells`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md)`(``fit``, ``expr``, cell_types ``=`` ``cts``,`\
`                            feature_scale ``=`` ``"simulated"``, min_coverage ``=`` ``0.9``)`\
`}`\
`projB`` ``<-`` ``project_cohort``(``fitA``, ``cohB``)`\
`projC`` ``<-`` ``project_cohort``(``fitA``, ``cohC``)`\
`projB``$``coverage``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"cell_type"``, ``"programs"``, ``"cells"``, ``"shared_genes"``, ``"coverage"``)``]`\
`#>   cell_type programs cells shared_genes coverage`\
`#> 1         T        2    24          400        1`\
`#> 2         B        1    24          400        1`

Coverage is complete (400/400 genes) — technical transportability holds
by construction here; on real data this table is where ineligible
type-cohorts are caught. Now the predeclared association on the
projected T-cell activities:

\
`assoc_projected`` ``<-`` ``function``(``proj``, ``meta``, ``ct`` ``=`` ``"T"``)`` ``{`\
`  ``A`` ``<-`` ``proj``$``activity`\
`  ``rows`` ``<-`` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``A``)``[`[`startsWith`](https://rdrr.io/r/base/startsWith.html)`(`[`rownames`](https://rdrr.io/r/base/colnames.html)`(``A``)``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``ct``, ``"_"``)``)``]`\
`  ``out`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``A``)``, ``function``(``p``)`` ``{`\
`    ``s`` ``<-`` ``A``[``rows``, ``p``]`\
`    ``if`` ``(`[`all`](https://rdrr.io/r/base/all.html)`(`[`is.na`](https://rdrr.io/r/base/NA.html)`(``s``)``)``)`` `[`return`](https://rdrr.io/r/base/function.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``NA_real_``, ``NA_real_``)``)`\
`    ``y`` ``<-`` ``meta``$``outcome``[`[`match`](https://rdrr.io/r/base/match.html)`(`[`sub`](https://rdrr.io/r/base/grep.html)`(``"^T_"``, ``""``, ``rows``)``, ``meta``$``observation_id``)``]`\
`    ``cf`` ``<-`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(`[`summary`](https://rdrr.io/r/base/summary.html)`(`[`lm`](https://rdrr.io/r/stats/lm.html)`(``y`` ``~`` ``s``)``)``)`\
`    `[`c`](https://rdrr.io/r/base/c.html)`(``cf``[``"s"``, ``1``]``, ``cf``[``"s"``, ``4``]``)`\
`  ``}``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``2``)``)`\
`  `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``program ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``A``)``, beta ``=`` ``out``[``1``, ``]``, p ``=`` ``out``[``2``, ``]``,`\
`             row.names ``=`` ``NULL``)`\
`}`\
`aB`` ``<-`` ``assoc_projected``(``projB``, ``cohB``$``meta``)`\
`aC`` ``<-`` ``assoc_projected``(``projC``, ``cohC``$``meta``)`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(`[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``cohort ``=`` ``"B"``, ``aB``)``, `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``cohort ``=`` ``"C"``, ``aC``)``)`\
`#>   cohort program       beta           p`\
`#> 1      B     T_1 -2.8684965 0.006144099`\
`#> 2      B     T_2  0.5826891 0.739876914`\
`#> 3      B     B_1         NA          NA`\
`#> 4      C     T_1 -3.1295870 0.003218153`\
`#> 5      C     T_2 -0.3099099 0.815471556`\
`#> 6      C     B_1         NA          NA`

One T-cell program shows the replicated outcome association in both
query cohorts (same sign as the discovery effect, p \< 0.05) — the
conserved P1 signal. Which discovery program is it, and is it the
planted one?

\
`sc`` ``<-`` `[`program_scores`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)`(``fitA``, format ``=`` ``"list"``)``$``T`\
`disc_assoc`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``sc``)``)``, ``function``(``k``)`` ``{`\
`  ``cf`` ``<-`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(`[`summary`](https://rdrr.io/r/base/summary.html)`(`[`lm`](https://rdrr.io/r/stats/lm.html)`(``cohA``$``meta``$``outcome`` ``~`` ``sc``[``, ``k``]``)``)``)``[``2``, `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``4``)``]`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``cf``[``1``]``, ``cf``[``2``]``)`\
`}``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``2``)``)`\
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``program ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sc``)``, beta ``=`` ``disc_assoc``[``1``, ``]``, p ``=`` ``disc_assoc``[``2``, ``]``)`\
`#>   program        beta            p`\
`#> 1     T_1 -0.77061789 2.472144e-05`\
`#> 2     T_2 -0.01467328 9.472151e-01`

The same discovery program carries the association in A and replicates
in B and C through the fixed basis — score-level transfer plus effect
replication, with the sign fixed by canonicalization within the
projected frame.

## Pathway 2: de-novo refit + matching

Cohort B is refit independently — same frozen feature panel, independent
factors — then aligned to the discovery programs by loading-space
matching:

\
`fitB`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``cohB``$``data``, max_factors ``=`` ``4``, seed ``=`` ``2``,`\
`                              preprocessing_spec ``=`` ``frozen``)`\
`fitB`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fitB``)`\
`m`` ``<-`` `[`match_programs`](https://mdmanurung.github.io/cellprograms/reference/match_programs.md)`(``fitA``, ``fitB``, ``"T"``, method ``=`` ``"correlation"``,`\
`                    space ``=`` ``"loadings"``, sim_threshold ``=`` ``0.5``)`\
`m``$``matches`\
`#>   program_a program_b    cosine`\
`#> 1       T_1       T_1 0.9641598`\
`#> 2       T_2       T_2 0.8973300`

The matched B program’s outcome association is computed on **B’s own
scores** (never on discovery coordinates):

\
`scB`` ``<-`` `[`program_scores`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)`(``fitB``, format ``=`` ``"list"``)``$``T`\
`hit`` ``<-`` ``m``$``matches``[``m``$``matches``$``cosine`` ``>`` ``0.5``, , drop ``=`` ``FALSE``]`\
`if`` ``(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``hit``)``)`` ``{`\
`  ``k`` ``<-`` `[`match`](https://rdrr.io/r/base/match.html)`(``hit``$``program_b``[``1``]``, `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``scB``)``)`\
`  `[`coef`](https://rdrr.io/r/stats/coef.html)`(`[`summary`](https://rdrr.io/r/base/summary.html)`(`[`lm`](https://rdrr.io/r/stats/lm.html)`(``cohB``$``meta``$``outcome`` ``~`` ``scB``[``, ``k``]``)``)``)``[``2``, ``]`\
`}`\
`#>     Estimate   Std. Error      t value     Pr(>|t|) `\
`#> -0.506379271  0.167712860 -3.019322857  0.006305575`

## Abstention accounting

Every program enters exactly one bucket per pathway; denominators are
always shown.

\
`disc_sign`` ``<-`` `[`sign`](https://rdrr.io/r/base/sign.html)`(``disc_assoc``[``1``, `[`which.min`](https://rdrr.io/r/base/which.min.html)`(``disc_assoc``[``2``, ``]``)``]``)``  ``# discovery effect direction`\
`accounting`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  pathway ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"fixed-basis (A->B)"``, ``"fixed-basis (A->C)"``, ``"de-novo (A vs B)"``)``,`\
`  discovery_programs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitA``$``scores``$``T``)``, `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitA``$``scores``$``T``)``, `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitA``$``scores``$``T``)``)``,`\
`  attempted ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(`[`startsWith`](https://rdrr.io/r/base/startsWith.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``projB``$``activity``)``, ``"T_"``)``)``,`\
`                `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`startsWith`](https://rdrr.io/r/base/startsWith.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``projC``$``activity``)``, ``"T_"``)``)``,`\
`                `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``m``$``matches``)``)``,`\
`  ineligible ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0L``, ``0L``, `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitA``$``scores``$``T``)`` ``-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``m``$``matches``)``)``,`\
`  replicated ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(``aB``$``p`` ``<`` ``0.05`` ``&`` `[`sign`](https://rdrr.io/r/base/sign.html)`(``aB``$``beta``)`` ``==`` ``disc_sign``, na.rm ``=`` ``TRUE``)``,`\
`                 `[`sum`](https://rdrr.io/r/base/sum.html)`(``aC``$``p`` ``<`` ``0.05`` ``&`` `[`sign`](https://rdrr.io/r/base/sign.html)`(``aC``$``beta``)`` ``==`` ``disc_sign``, na.rm ``=`` ``TRUE``)``,`\
`                 ``NA_integer_``)``,`\
`  stringsAsFactors ``=`` ``FALSE``)`\
`accounting`\
`#>              pathway discovery_programs attempted ineligible replicated`\
`#> 1 fixed-basis (A->B)                  2         2          0          1`\
`#> 2 fixed-basis (A->C)                  2         2          0          1`\
`#> 3   de-novo (A vs B)                  2         2          0         NA`

And the cohort-specific truth: P3 exists only in cohort C. When C is the
discovery cohort, its extra program must **abstain** in A and B — no
projection target, no match. That is the correct outcome, not a failure:

\
`fitC`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``cohC``$``data``, max_factors ``=`` ``5``, seed ``=`` ``3``,`\
`                              feature_scale ``=`` ``"simulated"``)`\
`fitC`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fitC``)`\
`m_CA`` ``<-`` `[`match_programs`](https://mdmanurung.github.io/cellprograms/reference/match_programs.md)`(``fitC``, ``fitA``, ``"T"``, method ``=`` ``"correlation"``,`\
`                       space ``=`` ``"loadings"``, sim_threshold ``=`` ``0.5``)`\
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``discovery_C_programs ``=`` `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitC``$``scores``$``T``)``,`\
`           matched_in_A ``=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``m_CA``$``matches``)``,`\
`           abstained ``=`` `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``fitC``$``scores``$``T``)`` ``-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``m_CA``$``matches``)``)`\
`#>   discovery_C_programs matched_in_A abstained`\
`#> 1                    3            2         1`

C discovers more T programs than A (it carries P3); the unmatched
remainder is the abstention count. The leave-one-cohort-out gate is:
conserved P1 replicates in every held-out cohort (shown above for B and
C with A as discovery; the rotations B-\>A,C and C-\>A,B use the same
code), and cohort-specific programs abstain rather than force a match.

## Negative controls

**(1) Permuted outcome.** Replication must collapse when the outcome is
permuted within cohort:

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``99``)`\
`perm_hits`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``200``)``, ``function``(``r``)`` ``{`\
`  ``meta_p`` ``<-`` ``cohB``$``meta`\
`  ``meta_p``$``outcome`` ``<-`` `[`sample`](https://rdrr.io/r/base/sample.html)`(``meta_p``$``outcome``)``   ``# shuffle values against fixed donor IDs`\
`  ``aB_perm`` ``<-`` ``assoc_projected``(``projB``, ``meta_p``)`\
`  `[`sum`](https://rdrr.io/r/base/sum.html)`(``aB_perm``$``p`` ``<`` ``0.05`` ``&`` `[`sign`](https://rdrr.io/r/base/sign.html)`(``aB_perm``$``beta``)`` ``==`` ``disc_sign``, na.rm ``=`` ``TRUE``)`\
`}``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)``)`\
[`mean`](https://rdrr.io/r/base/mean.html)`(``perm_hits``)``  ``# expectation under the null: ~0.05 per non-degenerate program`\
`#> [1] 0.07`

**(2) Scrambled gene IDs.** A query whose feature IDs do not match the
reference must abstain loudly, not silently project:

\
`bad`` ``<-`` ``cohB`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``bad``$``data``$``pseudobulk``$``T``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"X"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_genes``)``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``bad``$``data``$``pseudobulk``$``B``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"X"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_genes``)``)`\
[`tryCatch`](https://rdrr.io/r/base/conditions.html)`(``project_cohort``(``fitA``, ``bad``)``,`\
`         error ``=`` ``function``(``e``)`` `[`conditionMessage`](https://rdrr.io/r/base/conditions.html)`(``e``)``)`\
`` #> [1] "\033[1m\033[22mCell type \033[34m\"T\"\033[39m: none of the fit's 400 genes are present in `expression`." ``

**(3) Nuisance-only effect.** P2 tracks batch in cohort B by
construction. A technical association must not present as replicated
biology: testing the batch term in cohort C (where no batch effect was
planted) finds nothing.

\
`sB`` ``<-`` ``projB``$``activity``[`[`paste0`](https://rdrr.io/r/base/paste.html)`(``"T_"``, ``cohB``$``meta``$``observation_id``)``, ``]`\
`sC`` ``<-`` ``projC``$``activity``[`[`paste0`](https://rdrr.io/r/base/paste.html)`(``"T_"``, ``cohC``$``meta``$``observation_id``)``, ``]`\
`batch_assoc`` ``<-`` ``function``(``S``, ``meta``)`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``S``)``)``, ``function``(``k``)`` ``{`\
`  ``s`` ``<-`` ``S``[``, ``k``]`\
`  ``if`` ``(`[`all`](https://rdrr.io/r/base/all.html)`(`[`is.na`](https://rdrr.io/r/base/NA.html)`(``s``)``)``)`` `[`return`](https://rdrr.io/r/base/function.html)`(``NA_real_``)``   ``# off-cell-type program column`\
`  ``cf`` ``<-`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(`[`summary`](https://rdrr.io/r/base/summary.html)`(`[`lm`](https://rdrr.io/r/stats/lm.html)`(``s`` ``~`` ``meta``$``batch``)``)``)`\
`  ``if`` ``(``"meta$batchb2"`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``cf``)``)`` ``cf``[``"meta$batchb2"``, ``4``]`` ``else`` ``NA_real_`\
`}``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)``)`\
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``program ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sB``)``,`\
`           p_batch_B ``=`` ``batch_assoc``(``sB``, ``cohB``$``meta``)``,`\
`           p_batch_C ``=`` ``batch_assoc``(``sC``, ``cohC``$``meta``)``)`\
`#>   program    p_batch_B p_batch_C`\
`#> 1     T_1 8.860339e-01 0.5783604`\
`#> 2     T_2 7.579566e-09 0.4969985`\
`#> 3     B_1           NA        NA`

The batch-loaded program is significant in B and null in C — flagged as
cohort-specific technical structure, excluded from the replication
claim.

## What to report, and limitations

A replication report from this workflow should state: cohorts and
programs attempted/projected/matched/ineligible/nonreplicated (the
accounting table); the predeclared association and replication rule;
coverage and reconstruction per transfer (technical transportability);
and the negative-control outcomes. Limitations: projected activities are
query-set-dependent empirical-Bayes estimates (see the projection
documentation); sign and scale are comparable only within a pathway, not
across de-novo fits; and with 24 donors per cohort the power calculus
here is illustrative, not a recommendation.
