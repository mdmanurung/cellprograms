#!/bin/bash
#SBATCH --job-name=envfix2
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cellprograms/benchmarks/biomni_replication/results/logs/env_fix2.%j.log
export R_LIBS_USER=/nonexistent
export PATH=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin:$PATH
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript -e 'options(repos="https://cloud.r-project.org", Ncpus=4); install.packages(c("ebnm","flashier"), lib=.Library)' || exit 1
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript -e 'library(flashier); library(ebnm); library(class); cat("OK\n")'
