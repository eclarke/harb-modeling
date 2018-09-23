#!/usr/bin/env Rscript
# 
# Autoregressive binomial model, v2.0.14
# 2018-09-19
# 
# Data preparation and RStan execution script with options to specify
# taxa and data parameters


# Parse arguments ---------------------------------------------------------
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--genus", help="Genus to test", required=TRUE)
parser$add_argument(
  "--species", 
  help="Species to test (optional, if missing will aggregate by genus)")
parser$add_argument(
  "--min_subjects", default=2,
  help="Min. num. of subjects to have received an antibiotic for it to be included")
parser$add_argument(
  "--min_med_days", default=2,
  help="Minimum value of the median number of days an antibiotic was given across subjects")
parser$add_argument(
  "--specimen_type", 
  help="Specimen type to use ('Sputum', 'Oral Swab', 'Stool Swab')", required=TRUE)
parser$add_argument(
  "--min_reads_subj", default=0,
  help="Minimum number of reads present across a subject's samples to include that subject")
parser$add_argument("--iter", default=5000, help="stan iterations")
parser$add_argument("--max_treedepth", default=18, help="stan max_treedepth")
parser$add_argument("--adapt_delta", default=0.99, help="stan adapt_delta")
parser$add_argument("--cores", default=4, help="stan cores")
parser$add_argument("--chains", default=4, help="stan chains")
parser$add_argument("--model", help="Stan file to run", required=TRUE)
parser$add_argument("--data_gen", help="data generation function to use", required=TRUE)
parser$add_argument(
  "--output_fp", help="custom directory to store output files")

args <- parser$parse_args()

if (!is.null(args$output_fp) && dir.exists(args$output_fp)) {
  stop(sprintf("output directory '%s' already exists.", args$output_fp))
} 

# Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(phyloseq)
  library(tidyverse)
  library(stringr)
  library(magrittr)
  library(rstan)
  library(tidybayes)
})
source(here("shared_functions.R"))
source(here("betancourt_utils.R"))
source(here("data_generation.R"))

rstan_options(auto_write = TRUE)

if (!exists(args$data_gen, mode="function")) {
  stop(sprintf("specified --data_gen function '%s' not found", args$data_gen))
}

if (is.null(args$output_fp)) {
  args$output_fp <- sprintf(
    "_%s_%s_%s_%s",
    paste(c(args$genus, args$species), collapse="-"),
    stringr::str_replace(args$specimen_type, " ", ""),
    tools::file_path_sans_ext(basename(args$model)),
    args$data_gen)
  cat(sprintf("\nOutput directory: '%s'\n\n", args$output_fp))
}

if (dir.exists(args$output_fp)) {
  stop(sprintf("output directory '%s' already exists.", args$output_fp))
} 

# Load data ---------------------------------------------------------------
message("Loading data, please stand by...")
suppressMessages({
  d2 <- load_data(.genus=args$genus, .species=args$species)
})

dat <- do.call(
  args$data_gen, list(d2, args$specimen_type, args$min_subjects, args$min_med_days, args$min_reads_subj))

dir.create(args$output_fp)

saveRDS(dat, file.path(args$output_fp, "data.rds"))

cat(sprintf(
  "\n##\n## %d subjects, %d antibiotics, %d specimens; %1.2f%% zeros in read counts\n##\n\n", 
  dat$n_subjects, dat$n_abx, dat$n_specimens, sum(dat$d_abx$reads==0)/length(dat$d_abx$reads)*100))

m <- stan(
  file=here(args$model), data=dat, cores=args$cores, iter=args$iter, chains=args$chains, 
  control=list(max_treedepth=as.integer(args$max_treedepth), adapt_delta=as.double(args$adapt_delta)))

check_all_diagnostics(m)

saveRDS(m, file.path(args$output_fp, "fitted_model.rds"))
