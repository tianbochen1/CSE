# Daily returns, rolling comparisons, and refitted SV model checks.
CORES <- 8L
BOOTSTRAP_SAMPLES <- 499L
APPLICATIONS <- c("EURUSD", "SP500")
LOCAL_CSV <- list(EURUSD="",SP500="")
script_argument <- grep("^--file=",commandArgs(FALSE),value=TRUE)
SOURCE_DIR <- if(length(script_argument)) dirname(normalizePath(
  sub("^--file=","",script_argument[1L]),mustWork=TRUE)) else "R"
ROOT <- dirname(SOURCE_DIR)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
source(file.path(SOURCE_DIR,"cse_core.R"))
source(file.path(SOURCE_DIR,"cse_applications.R"))
CFG <- list(trunc_tol=1e-5,gh_order=512L,bank_nodes=81L,bank_max_nodes=641L,
  bank_tol=2e-6,quadrature_tol=2e-6,cpl_grid=49L,cpl_gh=32L,
  lower=c(.005,.01),upper=c(.995,3),nstarts=4L,maxit=300L)
DATA_CFG <- list(data_dir=file.path(ROOT,"data"),local_csv=LOCAL_CSV,refresh=FALSE,
  start=as.Date("2017-01-01"),end=as.Date("2025-12-31"),train_n=1024L,test_n=256L)
INPUT <- setNames(lapply(APPLICATIONS,cse_load_application,cfg=DATA_CFG),APPLICATIONS)
EVALUATION <- cse_eval_spec()
Q_MULTIPLIER <- 16L
F_MULTIPLIER <- 32L
levels <- unique(c(EVALUATION$tau,unlist(lapply(c(512L,1024L,
  vapply(INPUT,function(z)nrow(z$data),integer(1))),function(n)
    cse_application_spec(n,Q_MULTIPLIER,F_MULTIPLIER)$tau))))
BANK <- cse_build_bank(levels,CFG)
OUTPUT <- file.path(ROOT,"results","applications")
dir.create(OUTPUT,recursive=TRUE,showWarnings=FALSE)
FOLDS <- do.call(rbind,lapply(seq_len(4L),function(i) do.call(rbind,
  lapply(c(512L,1024L),function(n) data.frame(fold=i,train_n=n,
    train_start=1024L+(i-1L)*256L-n+1L,train_end=1024L+(i-1L)*256L,
    test_start=1025L+(i-1L)*256L,test_end=1024L+i*256L)))))
POINTS <- setNames(lapply(APPLICATIONS,function(a) {
  y<-INPUT[[a]]$data$y
  cse_fit_real(y,"CSE",cse_application_spec(length(y),Q_MULTIPLIER,F_MULTIPLIER),BANK,CFG)
}),APPLICATIONS)
for(a in APPLICATIONS) {
  write.csv(INPUT[[a]]$audit,file.path(OUTPUT,paste0(a,"_data_summary.csv")),row.names=FALSE)
  write.csv(cse_fit_row(POINTS[[a]]),file.path(OUTPUT,paste0(a,"_point_estimate.csv")),row.names=FALSE)
  y<-INPUT[[a]]$data$y
  rows<-list()
  for(i in seq_len(nrow(FOLDS))) {
    f<-FOLDS[i,];train<-y[f$train_start:f$train_end];test<-y[f$test_start:f$test_end]
    spec<-cse_application_spec(length(train),Q_MULTIPLIER,F_MULTIPLIER)
    pilot<-cse_fit_real(train,"CSE",spec,BANK,CFG)
    for(method in c("CSE","TD2","TD16","QML","CPL","Independence")) {
      fit<-if(method=="CSE")pilot else if(method=="Independence")NULL else
        cse_fit_real(train,method,spec,BANK,CFG,extra=pilot$theta)
      scores<-cse_validation_score(test,if(is.null(fit))NULL else fit$theta,EVALUATION,BANK)
      rows[[length(rows)+1L]]<-data.frame(fold=f$fold,train_n=f$train_n,method=method,
        phi=if(is.null(fit))NA_real_ else fit$theta[1],v=if(is.null(fit))NA_real_ else fit$theta[2],
        boundary=if(is.null(fit))FALSE else fit$boundary,as.list(scores))
    }
  }
  blocks<-do.call(rbind,rows)
  summary<-aggregate(cbind(spectral_real,spectral_complex,lag_discrepancy,boundary)~train_n+method,blocks,mean)
  write.csv(blocks,file.path(OUTPUT,paste0(a,"_window_blocks.csv")),row.names=FALSE)
  write.csv(summary,file.path(OUTPUT,paste0(a,"_window_summary.csv")),row.names=FALSE)
}
# Ties are assessed descriptively for the rounded exchange-rate observations.
y<-INPUT$EURUSD$data$y
spec<-cse_application_spec(length(y),Q_MULTIPLIER,F_MULTIPLIER)
base<-cse_fit_row(POINTS$EURUSD);base$rank_rule<-"upper";ties<-list(base)
set.seed(20261020L)
rank_function <- cse_ranks
for(i in 0:20) {
  rule<-if(i==0L)"average" else "random"
  z<-rank(y,ties.method=rule)
  cse_ranks <- function(z) z/length(z)
  fit<-cse_fit_real(z,"CSE",spec,BANK,CFG,extra=POINTS$EURUSD$theta)
  cse_ranks <- rank_function
  row<-cse_fit_row(fit);row$rank_rule<-if(i==0L)rule else paste0(rule,"_",i)
  ties[[length(ties)+1L]]<-row
}
write.csv(do.call(rbind,ties),file.path(OUTPUT,"EURUSD_tie_sensitivity.csv"),row.names=FALSE)
# The continuous SV bootstrap is used for the equity model checks.
RNGkind("L'Ecuyer-CMRG");set.seed(20261020L)
for(i in seq_len(BOOTSTRAP_SAMPLES)) .Random.seed<-parallel::nextRNGStream(.Random.seed)
SEEDS<-vector("list",BOOTSTRAP_SAMPLES)
for(i in seq_len(BOOTSTRAP_SAMPLES)) {SEEDS[[i]]<-.Random.seed;.Random.seed<-parallel::nextRNGStream(.Random.seed)}
cse_bootstrap_worker<-function(i) {
  assign(".Random.seed",SEEDS[[i]],envir=.GlobalEnv)
  y<-cse_generate(POINTS$SP500$theta,cse_innovations(nrow(INPUT$SP500$data)))
  fit<-cse_fit_real(y,"CSE",cse_application_spec(length(y),Q_MULTIPLIER,F_MULTIPLIER),BANK,CFG,
                    extra=POINTS$SP500$theta)
  stat<-cse_component_checks(y,fit$theta,EVALUATION,BANK)
  list(fit=cse_fit_row(fit),statistics=stat$statistics,empirical=stat$empirical)
}
if(CORES==1L) BOOT<-lapply(seq_len(BOOTSTRAP_SAMPLES),cse_bootstrap_worker) else
if(.Platform$OS.type=="windows") {
  cluster<-parallel::makePSOCKcluster(CORES)
  variables<-unique(c(grep("^cse_",ls(envir=.GlobalEnv),value=TRUE),
    "BANK","CFG","INPUT","POINTS","SEEDS","EVALUATION","Q_MULTIPLIER","F_MULTIPLIER"))
  parallel::clusterExport(cluster,variables,envir=.GlobalEnv)
  BOOT<-tryCatch(parallel::parLapplyLB(cluster,seq_len(BOOTSTRAP_SAMPLES),cse_bootstrap_worker),
                 finally=parallel::stopCluster(cluster))
} else BOOT<-parallel::mclapply(seq_len(BOOTSTRAP_SAMPLES),cse_bootstrap_worker,
                               mc.cores=CORES,mc.preschedule=FALSE)
stopifnot(!any(vapply(BOOT,inherits,logical(1),"try-error")))
observed<-cse_component_checks(INPUT$SP500$data$y,POINTS$SP500$theta,EVALUATION,BANK)
statistics<-do.call(rbind,lapply(BOOT,`[[`,"statistics"))
checks<-data.frame(check=names(observed$statistics),statistic=as.numeric(observed$statistics),
  bootstrap_p=vapply(names(observed$statistics),function(name)
    (1+sum(statistics[,name]>=observed$statistics[[name]]))/(BOOTSTRAP_SAMPLES+1),numeric(1)))
write.csv(checks,file.path(OUTPUT,"SP500_model_checks.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(BOOT,`[[`,"fit")),file.path(OUTPUT,"SP500_bootstrap_estimates.csv"),row.names=FALSE)
saveRDS(list(observed=observed,bootstrap=BOOT,settings=CFG),file.path(OUTPUT,"SP500_spectral_checks.rds"))
