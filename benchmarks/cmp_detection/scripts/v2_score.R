## v2_score.R -- PLAN v2 scoring: scorecard + breadth diagnostics + paired bootstrap vs reference, on common datasets.
## Rscript v2_score.R <out.csv> <seed_expr e.g. 1:15> <ref_label=ref.rds[,ref_nullx.rds]> <label=path.rds[,nullx.rds]> ...
suppressPackageStartupMessages(library(data.table))
source("/mnt/shared-workspace/cmp/bench/config.R")
POS_SCEN <- c("cmp_distinct", "cmp_mixed", "context", "ambient", "lr", "weak")
a <- commandArgs(TRUE); out <- a[1]; seeds <- eval(parse(text = a[2]))
spec <- lapply(a[-(1:2)], function(x) { p <- strsplit(x, "=")[[1]]; f <- strsplit(p[2], ",")[[1]]; list(label = p[1], dev = f[1], nullx = if (length(f) > 1) f[2] else NA) })
seed_raw <- function(s) s %% 100000   # truth$seed carries a base/scenario offset
ld <- function(sp) {
  e <- readRDS(sp$dev); s <- as.data.table(e$summary); g <- as.data.table(e$signals)
  s[, sr := seed_raw(seed)]; g[, sr := seed_raw(seed)]
  s <- s[sr %in% seeds]; g <- g[sr %in% seeds]
  nx <- if (!is.na(sp$nullx) && file.exists(sp$nullx)) { x <- as.data.table(readRDS(sp$nullx)$summary); x } else NULL
  list(s = s, g = g, nx = nx)
}
D <- lapply(spec, ld); names(D) <- vapply(spec, `[[`, "", "label")
## common positive datasets across all entries (paired comparisons must use identical datasets)
keyf <- function(g) unique(g[scenario %in% POS_SCEN & mode == "marginal", paste(base, scenario, sr)])
common <- Reduce(intersect, lapply(D, function(d) keyf(d$g)))
brd <- function(sig) fifelse(grepl("R20", sig), "R20", fifelse(grepl("R50", sig), "R50", fifelse(sig == "CTX", "CTX", "R100")))
cp <- function(x) if (!length(x)) c(NA, NA) else stats::binom.test(sum(x), length(x))$conf.int[1:2]
ref <- D[[1]]
rows <- lapply(names(D), function(nm) {
  d <- D[[nm]]
  g <- d$g[scenario %in% POS_SCEN & mode == "marginal" & paste(base, scenario, sr) %in% common]
  g[, br := brd(signal)]
  o <- list(label = nm, n_datasets = length(common), F1 = mean(g$f1), upstream_recall = mean(g$oracle))
  for (b in c("param", "semi")) {
    gb <- g[base == b]; sb <- d$s[base == b]
    o[[paste0(b, "_F1")]] <- mean(gb$f1); o[[paste0(b, "_recall")]] <- mean(gb$oracle)
    for (r in c("R20", "R50", "R100", "CTX")) { x <- gb[br == r]; o[[paste0(b, "_F1_", r)]] <- mean(x$f1); o[[paste0(b, "_miss_", r)]] <- mean(x$f1 == 0) }
    nul <- sb[scenario == "null", any_false]
    if (!is.null(d$nx)) nul <- c(nul, d$nx[base == b & scenario == "null", any_false])
    o[[paste0(b, "_null_n")]] <- length(nul); o[[paste0(b, "_null_fwer")]] <- mean(nul); o[[paste0(b, "_null_cp_hi")]] <- cp(nul)[2]
    o[[paste0(b, "_ambient_decoy")]] <- mean(sb[scenario == "ambient", any_decoy])
    o[[paste0(b, "_within_false")]] <- mean(sb[scenario == "context" & mode == "within", any_false])
    o[[paste0(b, "_genes_only_false")]] <- mean(sb[scenario == "genes_only", any_false])
  }
  ## paired bootstrap vs reference on per-dataset mean F1
  agg <- function(gg) gg[scenario %in% POS_SCEN & mode == "marginal" & paste(base, scenario, sr) %in% common, .(f1 = mean(f1)), by = .(base, scenario, sr)]
  m <- merge(agg(ref$g), agg(d$g), by = c("base", "scenario", "sr"))
  dd <- m$f1.y - m$f1.x; set.seed(1); bs <- replicate(2000, mean(sample(dd, replace = TRUE)))
  o$diff_vs_ref <- mean(dd); o$diff_lo <- unname(quantile(bs, .025)); o$diff_hi <- unname(quantile(bs, .975))
  o$diff_param <- mean(dd[m$base == "param"]); o$diff_semi <- mean(dd[m$base == "semi"])
  gmax <- max(unlist(o[grep("(null_fwer|ambient_decoy|within_false|genes_only_false)$", names(o))]), na.rm = TRUE)
  o$max_guardrail <- gmax
  as.data.table(lapply(o, function(v) if (is.numeric(v)) round(v, 4) else v))
})
res <- rbindlist(rows, fill = TRUE)
fwrite(res, out)
print(res[, .(label, n_datasets, F1, upstream_recall, diff_vs_ref, diff_lo, diff_hi, param_null_fwer, semi_null_fwer, max_guardrail,
              param_miss_R20, semi_miss_R20, param_miss_R100)])
