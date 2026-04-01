args <- commandArgs(trailingOnly = TRUE)

bulk <- args[1]
python_home <- args[2]
results_path <- args[3]

source("SCCAF_D.R")

results <- SCCAF_D(
    bulk = bulk,
    python_home = python_home
)

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(
    results,
    file = file.path(results_path, "sccaf-d_results.rds")
)