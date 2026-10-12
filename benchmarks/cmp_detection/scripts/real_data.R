## real_data.R -- fit per-cell-type programs on a real cohort and run the frozen CMP finalist + annotations.
## Rscript real_data.R <stephenson_pb|combat_pb> [n_perm=999]
source("/mnt/shared-workspace/cmp/bench/common.R")
for (f in c("coordination.R", "cmp-annotate.R")) sys.source(file.path(CMP_ROOT, "pkg/R", f), envir = globalenv())
a <- commandArgs(TRUE); name <- a[1]; n_perm <- if (length(a) > 1) as.integer(a[2]) else 999L
covs <- strsplit(Sys.getenv("CMP_COV", "institute"), ",")[[1]]; tag <- if (length(covs) > 1) paste0("_", paste(covs[-1], collapse = "_")) else ""
kmax <- as.integer(Sys.getenv("CMP_KMAX", "10")); if (kmax != 10L) tag <- paste0(tag, "_K", kmax)   # v2: program cap
out <- file.path("/mnt/results/tusoskill/cmp_v1/real", paste0(name, tag)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, f) { tmp <- tempfile(fileext = ".csv"); write.csv(x, tmp, row.names = FALSE); system2("cp", c(tmp, file.path(out, f))) }
wrds <- function(x, f) { tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp); system2("cp", c(tmp, f)) }

b <- readRDS(file.path(CMP_ROOT, "data", name, "semi_base.rds"))
Y <- lapply(b$cts, `[[`, "Y")
ids <- sort(unique(unlist(lapply(Y, rownames))))
dm <- b$donor_meta; dm <- dm[match(ids, dm$donor_id), ]
if (name == "stephenson_pb") {
  grp <- dm$Status_on_day_collection_summary
  sev <- c(Healthy = 0, Asymptomatic = 1, Mild = 2, Moderate = 3, Severe = 4, Critical = 5)[grp]
  grp2 <- ifelse(grepl("^LPS", grp), "LPS", ifelse(grp == "Non_covid", "Non_covid", ifelse(grp == "Healthy", "Healthy", "COVID")))
} else {
  grp <- dm$Source
  sev <- c(HV = 0, COVID_HCW_MILD = 1, COVID_MILD = 2, COVID_SEV = 3, COVID_CRIT = 4)[grp]
  grp2 <- ifelse(grepl("^COVID", grp), "COVID", grp)
}
meta <- data.frame(observation_id = ids, institute = factor(dm$site), group = grp2, group_fine = grp,
                   severity = unname(sev), sex = dm$sex, stringsAsFactors = FALSE)
cat("donors", length(ids), "\n"); print(table(meta$institute, meta$group))

ff <- file.path(CMP_ROOT, "real", paste0(name, tag, "_fit.rds")); dir.create(dirname(ff), showWarnings = FALSE)
if (file.exists(ff) && file.size(ff) > 0) fit <- readRDS(ff) else {
  x <- as_cell_program_data(Y, sample_metadata = meta)
  t0 <- Sys.time()
  fit <- canonicalize_programs(fit_celltype_programs(x, loading_prior = "point_laplace", var_type = 1, features = "variable",
                                                     n_variable_genes = 2000L, covariate_mode = "residualize", seed = 1L,
                                                     max_factors = kmax, covariates = covs))
  cat("fit time", format(Sys.time() - t0), "\n")
  wrds(fit, ff)
}
if (is.null(fit$pseudobulk)) fit$pseudobulk <- Y
cat("programs per ct:\n"); print(vapply(fit$scores, ncol, 1L))
prof <- lapply(b$cts, function(z) log2(z$cpm_mean + 1))

## frozen finalist: node test only, tau 2, link 3, alpha 0.05; ambient flagged (and dropped for comparison)
cf <- file.path(CMP_ROOT, "real", paste0(name, tag, "_cmp.rds"))
res <- if (file.exists(cf) && file.size(cf) > 0) readRDS(cf) else list(
  flag   = find_coordinated_programs(fit, strata = "institute", ambient = "flag", profiles = prof, n_perm = n_perm),
  drop   = find_coordinated_programs(fit, strata = "institute", ambient = "drop", profiles = prof, n_perm = n_perm),
  within = find_coordinated_programs(fit, strata = "institute", context = "group", mode = "within", ambient = "flag",
                                     profiles = prof, n_perm = n_perm))
for (k in names(res)) { cat("\n==", k, "==\n"); print(res[[k]]) }
if (!file.exists(cf)) wrds(res, cf)

hall <- read.csv(file.path(CMP_ROOT, "res/hallmark.csv")); sets <- split(hall$gene_symbol, hall$gs_name)
lr <- read.csv(file.path(CMP_ROOT, "res/lr_consensus.csv"))
for (k in names(res)) {
  r <- res[[k]]
  if (!nrow(r$cmps)) next
  cm <- merge(r$cmps, r$ambient, by = c("cell_type", "program_id"), all.x = TRUE)
  wcsv(cm[order(cm$cmp, cm$cell_type), ], sprintf("%s_cmps.csv", k))
  wcsv(r$edges, sprintf("%s_edges.csv", k))
  wcsv(annotate_cmp_themes(fit, r, sets), sprintf("%s_themes.csv", k))
  wcsv(cmp_gene_overlap(fit, r), sprintf("%s_gene_overlap.csv", k))
  wcsv(annotate_cmp_lr(fit, r, lr, n_perm = 999L), sprintf("%s_lr.csv", k))
  ca <- rbind(cbind(variable = "severity", cmp_context_association(fit, r, "severity", covariates = "institute")),
              cbind(variable = "group", cmp_context_association(fit, r, "group", covariates = "institute")))
  wcsv(ca, sprintf("%s_context_assoc.csv", k))
}
wcsv(res$flag$nodes, "flag_nodes.csv")

## Site replication (score level; programs learned on all donors): edge correlations within each large site
if (name == "stephenson_pb" && nrow(res$flag$edges)) {
  e <- res$flag$edges; st <- setNames(as.character(meta$institute), meta$observation_id)
  rep <- do.call(rbind, lapply(c("Ncl", "Cambridge"), function(s) {
    do.call(rbind, lapply(seq_len(nrow(e)), function(i) {
      Sa <- fit$scores[[e$ct_a[i]]][, e$prog_a[i]]; Sb <- fit$scores[[e$ct_b[i]]][, e$prog_b[i]]
      d <- intersect(names(Sa)[st[names(Sa)] == s], names(Sb))
      ct <- suppressWarnings(stats::cor.test(Sa[d], Sb[d], method = "spearman"))
      data.frame(site = s, cmp = e$cmp[i], ct_a = e$ct_a[i], prog_a = e$prog_a[i], ct_b = e$ct_b[i], prog_b = e$prog_b[i],
                 full_sign = e$sign[i], rho = unname(ct$estimate), p = ct$p.value, n = length(d))
    }))
  }))
  wcsv(rep, "site_replication_edges.csv")
  ## independent detection within each site
  det <- do.call(rbind, lapply(c("Ncl", "Cambridge"), function(s) {
    keep <- meta$observation_id[meta$institute == s]
    fs <- fit; fs$observation_ids <- keep; fs$sample_metadata <- fit$sample_metadata[fit$sample_metadata$observation_id %in% keep, ]
    fs$scores <- lapply(fit$scores, function(S) S[intersect(rownames(S), keep), , drop = FALSE])
    r <- find_coordinated_programs(fs, ambient = "flag", profiles = prof, n_perm = n_perm, min_shared = 20L, min_donors = 20L)
    if (!nrow(r$cmps)) return(NULL)
    cbind(site = s, r$cmps)
  }))
  if (!is.null(det)) wcsv(det, "site_detection.csv")
  cat("\nsite replication summary\n")
  print(aggregate(cbind(concordant = sign(rep$rho) == rep$full_sign, nominal = sign(rep$rho) == rep$full_sign & rep$p < 0.05) ~ site,
                  rep, mean))
}
cat("done\n")
