#!/bin/bash
#SBATCH --job-name=envfix
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cellprograms/benchmarks/biomni_replication/results/logs/env_fix.%j.log
export MAMBA_ROOT_PREFIX=/exports/para-lipg-hpc/mdmanurung/tools/micromamba/root
MM=/exports/para-lipg-hpc/mdmanurung/tools/micromamba/bin/micromamba
cd /exports/para-lipg-hpc/mdmanurung/cellprograms/benchmarks/biomni_replication
$MM create -y -n cellprograms-r -f env_fix.yml || exit 1
RS=$MAMBA_ROOT_PREFIX/envs/cellprograms-r/bin/Rscript
export PATH=$MAMBA_ROOT_PREFIX/envs/cellprograms-r/bin:$PATH
$RS -e 'options(repos="https://cloud.r-project.org", Ncpus=4); install.packages(c("ebnm","flashier"))' || exit 1
$RS -e 'library(flashier); library(ebnm); cat("OK\n")'
