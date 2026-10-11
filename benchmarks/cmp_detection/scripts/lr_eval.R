## lr_eval.R -- Tier-2 LR annotation benchmark with the package functions on the lr scenario.
## For each lr-scenario fit: call CMPs (frozen finalist, ambient flag), annotate LR, and score against planted pairs.
## Variants: raw (p < 0.05) and donor-checked (raw + donor Spearman of L in sender vs R in receiver > 0 with p_cor < 0.05).
## Rscript lr_eval.R <split dev|test> <cores> <outcsv>
source("/mnt/shared-workspace/cmp/bench/common.R"); source("/mnt/shared-workspace/cmp/bench/config.R")
for (f in c("coordination.R", "cmp-annotate.R")) sys.source(file.path(CMP_ROOT, "pkg/R", f), envir = globalenv())
suppressPackageStartupMessages(library(parallel))
a <- commandArgs(TRUE); split <- a[1]; cores <- as.integer(a[2]); outcsv <- a[3]; ctrl <- if (length(a) > 3) a[4] else "fwer"
seeds <- if (split == "dev") DEV_SEEDS else 1001:1050
lr <- read.csv(file.path(CMP_ROOT, "res/lr_consensus.csv"), stringsAsFactors = FALSE)
files <- unlist(lapply(c("param", "semi"), function(b) file.path(CMP_ROOT, "fits/main", b, "lr", sprintf("s%04d.rds", seeds))))
files <- files[file.exists(files)]
one <- function(f) {
  x <- readRDS(f); fit <- x$fit
  if (is.character(fit)) return(NULL)
  d <- gen_cmp(x$truth$base, "lr", x$truth$seed - ((match(x$truth$base, c("param", "semi")) - 1L) * 5e6 + match("lr", CMP_SCEN) * 1e5), x$amp)
  fit$pseudobulk <- d$x$pseudobulk %||% d$x[[1]]
  fit$sample_metadata$strata_ <- unname(x$strata[fit$sample_metadata$observation_id])
  r <- find_coordinated_programs(fit, strata = "strata_", ambient = "none", n_perm = 199)
  pl <- x$truth$lr; key_pl <- paste(pl$sender, pl$receiver, pl$ligand, pl$receptor)
  sig_genes <- unique(unlist(lapply(x$truth$signals, function(s) unlist(s$genes))))
  if (!nrow(r$cmps)) return(data.frame(file = f, base = x$truth$base, n_planted = nrow(pl), n_hit = 0, n_tp = 0, n_tp_dc = 0,
                                       n_hit_dc = 0, n_fp_ingenes = 0, n_fp_ingenes_dc = 0))
  h <- annotate_cmp_lr(fit, r, lr, n_perm = 499, control = ctrl, top = 1000L)
  if (!nrow(h)) h$donor_cor <- numeric()
  h$cs <- sub("\\|.*", "", h$sender); h$cr <- sub("\\|.*", "", h$receiver)
  ## donor-check p-value
  h$dc_p <- vapply(seq_len(nrow(h)), function(i) {
    Ys <- fit$pseudobulk[[h$cs[i]]]; Yr <- fit$pseudobulk[[h$cr[i]]]; dd <- intersect(rownames(Ys), rownames(Yr))
    if (!h$ligand[i] %in% colnames(Ys) || !h$receptor[i] %in% colnames(Yr)) return(NA_real_)
    suppressWarnings(stats::cor.test(Ys[dd, h$ligand[i]], Yr[dd, h$receptor[i]], method = "spearman")$p.value)
  }, 1)
  tp <- paste(h$cs, h$cr, h$ligand, h$receptor) %in% key_pl
  dc <- !is.na(h$donor_cor) & h$donor_cor > 0 & !is.na(h$dc_p) & h$dc_p < 0.05
  ing <- h$ligand %in% sig_genes & h$receptor %in% sig_genes
  data.frame(file = f, base = x$truth$base, n_planted = nrow(pl), n_hit = nrow(h), n_tp = sum(tp), n_tp_dc = sum(tp & dc),
             n_hit_dc = sum(dc), n_fp_ingenes = sum(!tp & ing), n_fp_ingenes_dc = sum(!tp & ing & dc))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
res <- do.call(rbind, mclapply(files, function(f) tryCatch(one(f), error = function(e) { message(f, ": ", conditionMessage(e)); NULL }),
                               mc.cores = cores, mc.preschedule = FALSE))
tmp <- tempfile(fileext = ".csv"); write.csv(res, tmp, row.names = FALSE); system2("cp", c(tmp, outcsv))
agg <- do.call(rbind, lapply(split(res, res$base), function(z) data.frame(
  base = z$base[1], n_files = nrow(z), recall = sum(z$n_tp) / sum(z$n_planted), recall_dc = sum(z$n_tp_dc) / sum(z$n_planted),
  precision = sum(z$n_tp) / max(1, sum(z$n_hit)), precision_dc = sum(z$n_tp_dc) / max(1, sum(z$n_hit_dc)),
  hits = sum(z$n_hit), hits_dc = sum(z$n_hit_dc), fp_in_program_genes = sum(z$n_fp_ingenes), fp_in_program_genes_dc = sum(z$n_fp_ingenes_dc))))
print(agg, row.names = FALSE)
