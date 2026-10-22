## common.R -- load repo package code + benchmark helpers on any machine
## CMP_ROOT env overrides the legacy layout; default = this repo checkout.
if (dir.exists("/workspace/.Rlib")) .libPaths(c("/workspace/.Rlib", .libPaths()))
HERE <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO_ROOT <- normalizePath(file.path(HERE, "..", ".."))
CMP_ROOT <- Sys.getenv("CMP_ROOT", unset = REPO_ROOT)
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
for (f in list.files(file.path(REPO_ROOT, "R"), pattern = "\\.R$", full.names = TRUE)) sys.source(f, envir = globalenv())
source(file.path(REPO_ROOT, "benchmarks/coordination/methods.R"))
source(file.path(HERE, "sim_cmp.R"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
