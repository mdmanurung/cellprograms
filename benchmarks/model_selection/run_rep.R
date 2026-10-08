## run_rep.R <base: param|semi> <scenario: S0|S1|S2|S3|S7> <rep> [amp] [outdir] [arms,comma]
## One dataset, every arm (paired). Writes <outdir>/rep<N>_ct.csv and rep<N>_pair.csv.
## Env MS_SMOKE=1 shrinks everything for a quick local check.
a <- commandArgs(TRUE); stopifnot(length(a) >= 3)
base <- a[1]; scenario <- a[2]; rep <- as.integer(a[3])
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_param.R", "sim_semi.R", "arms.R")) source(file.path(here, f))
stopifnot(base %in% c("param", "semi"), scenario %in% c("S0", "S1", "S2", "S3", "S7"))
if (nzchar(Sys.getenv("MS_SMOKE"))) options(ms.G = 600L, ms.nvar = 300L, ms.nperm = 49L, ms.nboot = 3L)

amp <- if (scenario == "S0") 0 else if (length(a) >= 4 && nzchar(a[4])) as.numeric(a[4]) else {
  tab <- utils::read.csv(file.path(here, "amplitudes.csv"), stringsAsFactors = FALSE)
  tab$amp[tab$base == base & tab$scenario == ifelse(scenario == "S7", "S1", scenario)]
}
stopifnot(length(amp) == 1, is.finite(amp))
outdir <- if (length(a) >= 5 && nzchar(a[5])) a[5] else file.path(here, "results", base, scenario)
arms <- if (length(a) >= 6 && nzchar(a[6])) strsplit(a[6], ",")[[1]] else ARMS
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- (base == "semi") * 5e6 + match(scenario, c("S0", "S1", "S2", "S3", "S7")) * 1e5 + rep
d <- if (base == "param") gen_param(scenario, amp, seed) else gen_semi(scenario, amp, seed, root)
res <- run_arms(d, arms)
tag <- data.frame(base = base, scenario = scenario, rep = rep, amp = amp, seed = seed)
write.csv(cbind(tag, res$ct), file.path(outdir, sprintf("rep%d_ct.csv", rep)), row.names = FALSE)
write.csv(cbind(tag, res$pair), file.path(outdir, sprintf("rep%d_pair.csv", rep)), row.names = FALSE)
