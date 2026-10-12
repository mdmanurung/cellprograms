## real_compare_K.R -- compare CMP calls between two program caps on the same cohort (same donors).
## Rscript real_compare_K.R <ref_run> <new_run> [<out_prefix>]
## Matching: |Pearson r| between CMP consensus donor scores (cmp_scores()); a new-run CMP is "retained" if its best
## |r| with any ref CMP is >= 0.5, otherwise "new". Ref CMPs with no new-run partner >= 0.5 are "lost".
## Per CMP: breadth, VDJ dominance (>= 50% of recurrent top genes match the V(D)J/Ig regex), ambient majority,
## severity association, and (Stephenson) site replication of its edges and independent re-detection per site.
suppressPackageStartupMessages(library(data.table))
source("/mnt/shared-workspace/cmp/pkg/R/coordination.R")
R <- "/mnt/shared-workspace/cmp/real"; O <- "/mnt/results/tusoskill/cmp_v1/real"
a <- commandArgs(TRUE); ref <- a[1]; new <- a[2]; pre <- if (length(a) >= 3) a[3] else paste0(ref, "_vs_", new)
VDJ <- "^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)"
## genes with common germline deletion / copy-number / HLA-haplotype expression variation (donor genotype markers)
GENO <- "^(GSTM1|GSTT1|GSTT2B?|UTS2|ERAP2|CCL3L1|CCL3L3|CCL4L2|HLA-DRB[345]|HLA-DQ[AB][12]|XIST|RPS4Y1|DDX3Y|UTY|KDM5D|EIF1AY)$"
IEG <- "^(FOS|FOSB|JUN|JUNB|EGR1|EGR2|EGR3|ATF3|IER2|IER3|NR4A1|DUSP1|PPP1R15A|KLF6)$"
ld <- function(n) {
  fit <- readRDS(file.path(R, paste0(n, "_fit.rds"))); cm <- readRDS(file.path(R, paste0(n, "_cmp.rds")))$flag$cmps
  s <- fread(file.path(O, n, "flag_cmp_summary.csv"))
  list(sc = cmp_scores(fit, cm), cm = as.data.table(cm), s = s, dir = file.path(O, n))
}
A <- ld(ref); B <- ld(new)
ids <- intersect(rownames(A$sc), rownames(B$sc))
r <- abs(suppressWarnings(cor(B$sc[ids, , drop = FALSE], A$sc[ids, , drop = FALSE], use = "pairwise.complete.obs")))
r[!is.finite(r)] <- 0
site_tab <- function(X) {
  f <- file.path(X$dir, "site_replication_edges.csv"); g <- file.path(X$dir, "site_detection.csv")
  if (!file.exists(f)) return(NULL)
  e <- fread(f)[, .(conc = mean(sign(rho) == full_sign), nom = mean(sign(rho) == full_sign & p < 0.05)), by = .(site, cmp)]
  e <- dcast(e, cmp ~ site, value.var = c("conc", "nom"))
  if (file.exists(g)) {
    d <- fread(g); full <- X$cm[, .(key = paste(cell_type, program_id)), by = cmp]
    det <- rbindlist(lapply(unique(full$cmp), function(k) {
      fk <- full[cmp == k]$key
      rbindlist(lapply(c("Ncl", "Cambridge"), function(s) {
        ds <- d[site == s]; if (!nrow(ds)) return(data.table(cmp = k, site = s, shared = 0L))
        sh <- ds[, .(shared = sum(paste(cell_type, program_id) %in% fk)), by = .(scmp = cmp)]
        data.table(cmp = k, site = s, shared = max(sh$shared))
      }))
    }))
    det <- dcast(det, cmp ~ site, value.var = "shared"); setnames(det, c("Ncl", "Cambridge"), c("det_members_Ncl", "det_members_Cambridge"))
    e <- merge(e, det, by = "cmp", all = TRUE)
  }
  e
}
ann <- function(X) {
  s <- copy(X$s); g <- strsplit(s$recurrent_top_genes, ",")
  s[, vdj_frac := sapply(g, function(z) mean(grepl(VDJ, trimws(z))))]
  s[, `:=`(vdj_dominated = vdj_frac >= 0.5, ambient_majority = n_ambient_flag > n_members / 2)]
  s[, n_geno := sapply(g, function(z) sum(grepl(GENO, trimws(z))))]
  s[, n_ieg := sapply(g, function(z) sum(grepl(IEG, trimws(z))))]
  st <- site_tab(X); if (!is.null(st)) s <- merge(s, st, by = "cmp", all.x = TRUE)
  s
}
SB <- ann(B); SA <- ann(A)
SB[, `:=`(best_ref = colnames(r)[apply(r, 1, which.max)][match(paste0("CMP", cmp), rownames(r))],
          best_ref_abs_r = round(apply(r, 1, max)[match(paste0("CMP", cmp), rownames(r))], 3))]
SB[, status := fifelse(best_ref_abs_r >= 0.5, "retained", "new")]
SA[, `:=`(best_new = rownames(r)[apply(r, 2, which.max)][match(paste0("CMP", cmp), colnames(r))],
          best_new_abs_r = round(apply(r, 2, max)[match(paste0("CMP", cmp), colnames(r))], 3))]
SA[, status := fifelse(best_new_abs_r >= 0.5, "retained", "lost")]
keep <- function(S) S[, intersect(c("cmp", "status", "best_ref", "best_ref_abs_r", "best_new", "best_new_abs_r", "n_members", "n_cell_types",
                                    "n_ambient_flag", "ambient_majority", "vdj_frac", "vdj_dominated", "n_geno", "n_ieg", "severity_est", "severity_p",
                                    "conc_Ncl", "conc_Cambridge", "nom_Ncl", "nom_Cambridge", "det_members_Ncl", "det_members_Cambridge",
                                    "top_themes", "recurrent_top_genes", "members"), names(S)), with = FALSE]
wr <- function(d, f) { t <- tempfile(fileext = ".csv"); fwrite(d, t); unlink(file.path(O, f)); system2("cp", c(t, file.path(O, f))) }
wr(keep(SB), paste0(pre, "_new_run_cmps.csv")); wr(keep(SA), paste0(pre, "_ref_run_cmps.csv"))
summ <- function(S, lab) data.table(run = lab, n_cmp = nrow(S), n_ge5ct = sum(S$n_cell_types >= 5), n_2to4ct = sum(S$n_cell_types < 5),
                                    vdj_dominated = sum(S$vdj_dominated), ambient_majority = sum(S$ambient_majority),
                                    sev_p_lt_1e3 = sum(S$severity_p < 1e-3, na.rm = TRUE), geno_marked = sum(S$n_geno >= 1), ieg_marked = sum(S$n_ieg >= 3))
out <- rbind(summ(SA, ref)[, `:=`(n_new = NA_integer_, n_new_narrow = NA_integer_, n_new_vdj = NA_integer_, n_new_amb = NA_integer_, n_new_geno = NA_integer_, n_new_sev = NA_integer_, n_lost = NA_integer_)],
             summ(SB, new)[, `:=`(n_new = sum(SB$status == "new"), n_new_narrow = sum(SB$status == "new" & SB$n_cell_types < 5),
                                  n_new_vdj = sum(SB$status == "new" & SB$vdj_dominated), n_new_amb = sum(SB$status == "new" & SB$ambient_majority), n_new_geno = sum(SB$status == "new" & SB$n_geno >= 1), n_new_sev = sum(SB$status == "new" & SB$severity_p < 0.05, na.rm = TRUE),
                                  n_lost = sum(SA$status == "lost"))])
wr(out, paste0(pre, "_summary.csv"))
print(out)
cat("\nnew CMPs:\n"); print(keep(SB)[status == "new", intersect(c("cmp", "n_cell_types", "best_ref_abs_r", "vdj_dominated", "ambient_majority", "severity_p",
                                                                  "conc_Ncl", "conc_Cambridge", "nom_Ncl", "nom_Cambridge", "det_members_Ncl", "det_members_Cambridge", "top_themes"), names(SB)), with = FALSE])
cat("\nlost ref CMPs:\n"); print(keep(SA)[status == "lost", .(cmp, n_cell_types, best_new_abs_r, vdj_dominated, ambient_majority, severity_p, top_themes)])
