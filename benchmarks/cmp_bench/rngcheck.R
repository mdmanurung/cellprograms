source("/mnt/shared-workspace/cmp/bench/common.R")
library(parallel)
f <- function(base, amp) mclapply(1:3, function(s) { d <- gen_cmp(base, "lr", s, amp); paste(d$truth$lr$ligand[1], paste(d$truth$signals$R50$members, collapse=","), round(d$x$pseudobulk[[1]][1,1],4)) }, mc.cores = 3)
print(unlist(f("param", 0.45))); print(unlist(f("semi", 0.875)))
