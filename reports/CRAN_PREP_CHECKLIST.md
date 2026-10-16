# CRAN preparation checklist (cellprograms)

Date: 2026-10-16. Submission is inherently manual (accounted in the plan); this
checklist is everything that can be verified or prepared ahead.

## Local state (verified this session)

- [x] R CMD build succeeds; vignettes build (sandbox has no pdflatex — final
      `--as-cran` must run on a machine with TeX; local checks use
      `--no-manual --ignore-vignettes` on the built tarball).
- [x] R CMD check on built tarball: 0 errors / 0 warnings / 0 notes
      (`_R_CHECK_FORCE_SUGGESTS_=false`; optional Suggests covered in CI).
- [x] Full testthat suite green (558 pass / 0 fail at v0.2.0).
- [x] Every export has a roxygen page with a runnable example and at least one
      test (repo law).
- [x] NEWS.md current through v0.2.0 (Stage B section added).
- [x] DESCRIPTION: Imports complete (digest added 2026-10-16); Version 0.2.0.

## Before submission (manual or next session)

- [ ] `R CMD check --as-cran` on a machine with pdflatex (clears the manual/
      vignette checks the sandbox cannot run).
- [ ] Win-builder + mac-builder pre-submission runs (CRAN's web forms).
- [ ] Title/description fields in DESCRIPTION pass CRAN's style scrutiny
      (no quotes, no "R" in title, sentence case).
- [x] `utils::citEntry`/CITATION file: inst/CITATION added 2026-10-16
      (verified via readCitationFile); CITATION.cff already present.
- [ ] License check: dependencies' licenses compatible (flashier, ebnm,
      flashier.utils, cli, digest — all MIT/GPL-compatible with the current
      LICENSE; confirm flashier.utils's license at release time).
- [x] URLs in docs resolve: `urlchecker::url_check()` run 2026-10-16 —
      all 6 URLs correct, 0 fixes needed.
- [ ] Check `R CMD check --as-cran` NOTE about "New submission" — expected.
- [x] Package-landing pages: pkgdown site rebuilt 2026-10-16 (commit
      7a3269a) — current with all 6 articles (incl. immune_aging_cross_cohort)
      and the full 46-export reference; live pages return 200.
- [ ] Decide: submit v0.2.0 as-is or after the batch-benchmark outcome (if the
      registered rule fires and R/ changes, bump to 0.2.1/0.3.0 first — never
      submit a version whose NEWS/check state is in flux).

## Submission steps (manual)

1. `R CMD build .` on the release commit; upload the tarball at
   https://cran.r-project.org/submit.html.
2. Reply to the maintainer confirmation email within the window.
3. Archive the submission receipt in this file with the date.
