# Estimation, identification, and feature comparisons.
REPLICATIONS <- 500L
CORES <- 8L
EXPERIMENTS <- c("E01", "E02", "E03", "E04", "E05", "E10", "E11")
script_argument <- grep("^--file=", commandArgs(FALSE), value=TRUE)
SOURCE_DIR <- if(length(script_argument)) dirname(normalizePath(
  sub("^--file=", "", script_argument[1L]),mustWork=TRUE)) else "R"
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
source(file.path(SOURCE_DIR,"cse_core.R"))
source(file.path(SOURCE_DIR,"cse_experiments.R"))
CFG <- cse_config(REPLICATIONS,CORES,EXPERIMENTS)
CFG$output <- file.path(dirname(SOURCE_DIR),"results","simulation")
cse_run(CFG,SOURCE_DIR)
