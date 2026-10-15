source("/mnt/shared-workspace/cmp/bench/common.R")
for (base in c("param","semi")) for (sc in CMP_SCEN) {
  t0 <- Sys.time()
  r <- tryCatch(gen_cmp(base, sc, 1, amp = 0.5), error = function(e) { cat("ERR", base, sc, conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  sg <- r$truth$signals
  cat(sprintf("%-5s %-12s cts %d donors %d signals [%s] lr %d  %.1fs\n", base, sc, length(r$x$expression %||% r$x[[1]]),
      length(r$strata), paste(vapply(sg, function(s) paste0(s$id, ":", s$type, ":", length(s$members)), ""), collapse=" "),
      nrow(r$truth$lr), as.numeric(Sys.time()-t0, units="secs")))
  if (sc == "lr") print(r$truth$lr)
}
str(r$x, max.level = 1)
