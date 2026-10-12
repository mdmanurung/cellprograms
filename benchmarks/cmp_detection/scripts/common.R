## common.R -- load repo package code + benchmark helpers on any machine
if (dir.exists("/workspace/.Rlib")) .libPaths(c("/workspace/.Rlib", .libPaths()))
CMP_ROOT <- "/mnt/shared-workspace/cmp"
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
for (f in list.files(file.path(CMP_ROOT, "repo/R"), pattern = "\\.R$", full.names = TRUE)) sys.source(f, envir = globalenv())
source(file.path(CMP_ROOT, "repo/benchmarks/coordination/methods.R"))
source(file.path(CMP_ROOT, "bench/sim_cmp.R"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
