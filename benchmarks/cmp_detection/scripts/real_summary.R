## real_summary.R -- one-row-per-CMP summary of real-data results. Rscript real_summary.R <name> <variant flag|drop|within>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); name <- a[1]; v <- a[2]
R <- "/mnt/shared-workspace/cmp"; out <- file.path("/mnt/results/tusoskill/cmp_v1/real", name)
fit <- readRDS(file.path(R, "real", paste0(name, "_fit.rds"))); res <- readRDS(file.path(R, "real", paste0(name, "_cmp.rds")))[[v]]
th <- fread(file.path(out, sprintf("%s_themes.csv", v))); ca <- fread(file.path(out, sprintf("%s_context_assoc.csv", v)))
go <- fread(file.path(out, sprintf("%s_gene_overlap.csv", v)))
lrf <- file.path(out, sprintf("%s_lr.csv", v)); lr <- if (file.exists(lrf) && file.size(lrf) > 5) fread(lrf) else data.table()
cm <- as.data.table(res$cmps)
sexg <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TSIX", "ZFY", "USP9Y", "TXLNGY")
techg <- c("HBB", "HBA1", "HBA2", "PPBP", "PF4", "MALAT1", "NEAT1", "FOS", "JUN", "DUSP1", "HSPA1A", "HSPA1B")
summ <- cm[, {
  tops <- lapply(seq_len(.N), function(i) {
    w <- fit$loadings[[cell_type[i]]][, program_id[i]] * sign[i]; names(sort(-w))[1:15] })
  allg <- unlist(tops); tab <- sort(table(allg), decreasing = TRUE)
  t1 <- th[cmp == .BY$cmp & !is.na(set)][order(-t)][, .N, by = set][order(-N)]
  .(n_members = .N, n_cell_types = uniqueN(cell_type), n_ambient_flag = sum(ambient_flag %in% TRUE),
    members = paste0(program_id, ifelse(sign < 0, "(-)", ""), collapse = ", "),
    recurrent_top_genes = paste(head(names(tab), 12), collapse = ","),
    sex_genes_in_top15 = sum(allg %in% sexg), tech_genes_in_top15 = sum(allg %in% techg),
    top_themes = paste(head(sprintf("%s(%d)", sub("HALLMARK_", "", t1$set), t1$N), 3), collapse = "; "))
}, by = cmp]
summ <- merge(summ, go[, .(cmp, overlap_class = class)], by = "cmp", all.x = TRUE)
sv <- ca[variable == "severity", .(cmp = as.integer(cmp), severity_est = round(estimate, 3), severity_p = signif(p, 3))]
summ <- merge(summ, sv, by = "cmp", all.x = TRUE)
gp <- ca[variable == "group", .(group_terms = paste(sprintf("%s:%.2f(p=%.2g)", sub("^x", "", term), estimate, p), collapse = "; ")), by = .(cmp = as.integer(cmp))]
summ <- merge(summ, gp, by = "cmp", all.x = TRUE)
if (nrow(lr)) summ <- merge(summ, lr[, .(n_lr = .N, lr_top = paste(head(sprintf("%s->%s %s-%s", sub("\\|.*", "", sender), sub("\\|.*", "", receiver), ligand, receptor), 3), collapse = "; ")), by = cmp], by = "cmp", all.x = TRUE)
fwrite(summ, file.path(tempdir(), "s.csv")); system2("cp", c(file.path(tempdir(), "s.csv"), file.path(out, sprintf("%s_cmp_summary.csv", v))))
print(summ[, .(cmp, n_cell_types, n_ambient_flag, recurrent_top_genes, sex_genes_in_top15, top_themes, overlap_class, severity_est, severity_p)])
