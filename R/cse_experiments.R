cse_config <- function(reps=500L, cores=8L,
                       experiments=c("E01","E02","E03","E04","E05","E10","E11")) {
  stopifnot(reps >= 1L, cores >= 1L)
  list(cores=as.integer(cores),experiments=experiments,
       reps=as.integer(reps),reps_by_experiment=list(),
       output="results/simulation",resume=TRUE,seed=20261018L,
       lower=c(.005,.01),upper=c(.98,3),nstarts=4L,
       maxit=200L,sim_maxit=100L,sim_parameter_tol=2e-4,
       trunc_tol=1e-5,gh_order=512L,bank_nodes=81L,bank_max_nodes=641L,
       bank_tol=2e-6,quadrature_tol=2e-6,
       covariance_paths=200L,weight_ridge=.02,
       cpl_grid=49L,cpl_gh=32L,make_plots=TRUE,checkpoint_every=32L)
}

cse_cases <- function(experiment,cfg) {
  main<-expand.grid(n=c(256L,512L,1024L,2048L),phi=c(.3,.6,.9),v=c(.25,1,2))
  grid<-function(nn,theta) {
    do.call(rbind,lapply(nn,function(n)data.frame(n=n,phi=theta[,1L],v=theta[,2L])))
  }
  regular<-rbind(c(.6,1),c(.9,1))
  ans<-switch(experiment,
    E01=grid(1024L,regular),
    E02=main,
    E03=grid(c(512L,1024L),regular),
    E04=grid(c(256L,512L,1024L),rbind(regular,c(.9,2))),
    E05=grid(c(256L,512L,1024L),regular),
    E06=grid(c(512L,1024L,2048L),rbind(c(.02,1),c(.05,1),c(.1,1),c(.7,.02),c(.7,.05),c(.7,.1))),
    E07=grid(c(1024L,2048L,4096L),rbind(c(.8,1.5),c(.9,2))),
    E08={x<-grid(c(512L,1024L),regular);do.call(rbind,lapply(c(0,.005,.01,.03),function(p)transform(x,contamination=p)))},
    E09=grid(c(256L,512L,1024L),regular),
    E10=grid(c(512L,1024L),regular),
    E11=grid(c(256L,512L,1024L),regular),
    E12=grid(c(512L,1024L),regular),stop("Unknown experiment."))
  if(!"contamination" %in% names(ans)) ans$contamination<-0
  ans$innovation<-if(experiment=="E12") "t5" else "normal"
  ans$experiment<-experiment
  ans$cell<-seq_len(nrow(ans))
  rownames(ans)<-NULL
  ans
}

cse_all_tau <- function(cases) {
  tau<-c(.1,.2,.35,.5,.9)
  for(n in unique(cases$n)) for(k in c("full","central","tail","diagonal","fixed"))
    tau<-c(tau,cse_spec(n,k)$tau)
  for(n in unique(cases$n)) tau<-c(tau,cse_spec(n,"full",2L)$tau)
  sort(unique(round(tau,12L)))
}

cse_task_seeds <- function(cases,cfg,reps) {
  RNGkind("L'Ecuyer-CMRG")
  tasks<-list()
  for(i in seq_len(nrow(cases))) {
    cc<-cases[i,,drop=FALSE]
    set.seed(cfg$seed+10000L*as.integer(sub("E","",cc$experiment))+cc$cell)
    state<-.Random.seed
    for(r in seq_len(reps)) {
      tasks[[length(tasks)+1L]]<-list(id=sprintf("%s_C%03d_R%05d",cc$experiment,cc$cell,r),
                                    case=cc,replicate=r,seed=state)
      state<-parallel::nextRNGStream(state)
    }
  }
  tasks
}

cse_fit_record <- function(method,fitfun,truth,task,bank,evalspec,margin="identity",setting="main",model_correct=TRUE) {
  start<-proc.time()[3L]
  error<-""
  warnings<-character()
  fit<-tryCatch(withCallingHandlers(fitfun(),warning=function(w){
    warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e){error<<-conditionMessage(e);NULL})
  elapsed<-unname(proc.time()[3L]-start)
  base<-data.frame(task=task$id,experiment=task$case$experiment,cell=task$case$cell,
                   replicate=task$replicate,n=task$case$n,phi=truth[1L],v=truth[2L],
                   contamination=task$case$contamination,innovation=task$case$innovation,
                   margin=margin,setting=setting,method=method,phi_hat=NA_real_,v_hat=NA_real_,
                   target_interpretation=if(model_correct)"copula_parameter"else"clean_parameter_sensitivity",
                   spectral_loss=NA_real_,elapsed=elapsed,fit_seconds=elapsed,pilot_seconds=0,
                   calls=NA_integer_,objective=NA_real_,
                   convergence=99L,boundary=NA,valid=FALSE,reused=FALSE,
                   error=error,warnings=paste(unique(warnings),collapse=" | "),stringsAsFactors=FALSE)
  if(!is.null(fit) && length(fit$theta)==2L && all(is.finite(fit$theta))) {
    base$phi_hat<-fit$theta[1L];base$v_hat<-fit$theta[2L]
    base$calls<-fit$calls;base$objective<-fit$objective
    base$convergence<-fit$convergence;base$boundary<-fit$boundary;base$valid<-TRUE
    if(model_correct) base$spectral_loss<-sum(evalspec$weights*(cse_model(fit$theta,evalspec,bank)-
                                                                cse_model(truth,evalspec,bank))^2)
  }
  base
}

cse_run_task <- function(task) {
  cfg<-CSE_CONTEXT$cfg;bank<-CSE_CONTEXT$bank
  assign(".Random.seed",task$seed,envir=.GlobalEnv)
  cc<-task$case;n<-cc$n;truth<-c(cc$phi,cc$v);experiment<-cc$experiment
  y<-cse_generate(truth,cse_innovations(n,cc$innovation))
  contamination_count<-0L
  if(cc$contamination>0) {
    affected<-runif(n)<cc$contamination
    contamination_count<-sum(affected)
    y[affected]<-y[affected]+30*sample(c(-1,1),sum(affected),replace=TRUE)
  }
  spec<-cse_spec(n)
  evalspec<-cse_spec(n,"full",2L)
  rows<-list();diag<-list(contamination_count=contamination_count)
  correct<-cc$contamination==0 && cc$innovation=="normal"
  add<-function(name,fun,margin="identity",setting="main",uses_pilot=FALSE) {
    row<-cse_fit_record(name,fun,truth,task,bank,evalspec,margin,setting,correct)
    if(uses_pilot) {
      row$pilot_seconds<-rows[[1L]]$elapsed
      row$elapsed<-row$elapsed+row$pilot_seconds
    }
    rows[[length(rows)+1L]]<<-row
    c(row$phi_hat,row$v_hat)
  }
  base_theta<-add("CSE",function()cse_fit_spectral(y,spec,bank,cfg))
  if(any(!is.finite(base_theta))) base_theta<-c(.5,1)
  basic<-function(yy=y,margin="identity") {
    add("TD2",function()cse_fit_td(yy,spec,c(1L,2L),bank,cfg),margin)
    add("TD16",function()cse_fit_td(yy,spec,c(1L,2L,4L,8L,16L),bank,cfg),margin)
    add("LSV_QML",function()cse_fit_qml(yy,cfg),margin)
  }
  if(experiment %in% c("E02","E06","E08","E11","E12")) basic()
  if(experiment=="E03") {
    basic()
    rank_rows<-do.call(rbind,rows)[c(1L,2L,3L),,drop=FALSE]
    for(margin in c("asinh","cubic","skew","cauchy")) {
      yy<-cse_transform(y,margin,truth[2L],bank$gh)
      equal<-all(is.finite(yy)) && identical(cse_ranks(yy),cse_ranks(y))
      diag[[paste0("invariance_",margin)]]<-list(exact_ranks=equal,
        max_feature_difference=if(equal) 0 else max(abs(cse_empirical(yy,spec)-cse_empirical(y,spec))))
      if(!equal) stop("A marginal transformation failed the finite-precision rank check.")
      copies<-rank_rows;copies$margin<-margin;copies$reused<-TRUE;copies$elapsed<-NA_real_;copies$calls<-NA_integer_
      copies$fit_seconds<-NA_real_;copies$pilot_seconds<-NA_real_
      for(k in seq_len(nrow(copies))) rows[[length(rows)+1L]]<-copies[k,,drop=FALSE]
      add("LSV_QML",function()cse_fit_qml(yy,cfg),margin)
    }
  }
  if(experiment=="E04") {
    for(kind in c("central","tail","diagonal")) {
      sp<-cse_spec(n,kind)
      add(paste0("CSE_",kind),function()cse_fit_spectral(y,sp,bank,cfg))
      D<-cse_jacobian(truth,sp,bank)
      diag[[kind]]<-svd(D*sqrt(sp$weights),nu=0,nv=0)$d
    }
  }
  if(experiment=="E05") for(kind in c("low","middle","high")) {
    sp<-cse_spec(n,kind)
    add(paste0("CSE_",kind),function()cse_fit_spectral(y,sp,bank,cfg))
    D<-cse_jacobian(truth,sp,bank)
    diag[[kind]]<-svd(D*sqrt(sp$weights),nu=0,nv=0)$d
  }
  if(experiment=="E07") {
    assign(".Random.seed",parallel::nextRNGSubStream(task$seed),envir=.GlobalEnv)
    ratios<-c(1L,5L,10L,20L)
    inn<-cse_innovations(n*max(ratios))
    for(ratio in ratios) {
      small<-list(initial=inn$initial,eta=head(inn$eta,n*ratio),eps=head(inn$eps,n*ratio))
      add(paste0("CSE_SIM",ratio),function()cse_fit_simulated(y,spec,cfg,small,base_theta),uses_pilot=TRUE)
    }
    diag$simulation_ratio<-ratios
  }
  if(experiment=="E09") {
    sp<-cse_spec(n,"fixed")
    add("CSE_fixed_I",function()cse_fit_spectral(y,sp,bank,cfg,W=diag(length(sp$weights))))
    assign(".Random.seed",parallel::nextRNGSubStream(task$seed),envir=.GlobalEnv)
    begin<-proc.time()[3L]
    weighting<-cse_covariance_weight(base_theta,n,sp,cfg$covariance_paths,cfg)
    setup_time<-unname(proc.time()[3L]-begin)
    add("CSE_fixed_W",function()cse_fit_spectral(y,sp,bank,cfg,W=weighting$W,extra=base_theta),uses_pilot=TRUE)
    rows[[length(rows)]]$elapsed<-rows[[length(rows)]]$elapsed+setup_time
    rows[[length(rows)]]$fit_seconds<-rows[[length(rows)]]$fit_seconds+setup_time
    diag$weighting<-weighting[c("eigenvalues","ridge")]
  }
  if(experiment=="E10") {
    sp<-cse_spec(n,"full",2L)
    add("CSE_fine",function()cse_fit_spectral(y,sp,bank,cfg))
    coarse<-spec;coarse$lambda<-pi*(seq_len(8L)-.5)/8
    coarse$weights<-rep(ifelse(coarse$pairs[,1L]==coarse$pairs[,2L],1,2),each=8L)
    coarse$weights<-coarse$weights/sum(coarse$weights)
    add("CSE_8freq",function()cse_fit_spectral(y,coarse,bank,cfg))
    diag$numerical_audit<-bank$audit
  }
  if(experiment=="E11") {
    add("CPL_Q",function()cse_fit_cpl(y,bank,cfg,extra=base_theta),uses_pilot=TRUE)
    add("CPL_Q_fine",function()cse_fit_cpl(y,bank,cfg,extra=base_theta,refined=TRUE),uses_pilot=TRUE)
  }
  list(id=task$id,seed=task$seed,status="completed",error="",
       estimates=do.call(rbind,rows),diagnostics=diag)
}

cse_safe_task <- function(task) {
  tryCatch(cse_run_task(task),error=function(e)
    list(id=task$id,seed=task$seed,status="failed",error=conditionMessage(e),
         case=task$case,replicate=task$replicate,estimates=NULL,diagnostics=NULL))
}

cse_summary <- function(d) {
  if(is.null(d)||!nrow(d)) return(data.frame())
  groups<-c("experiment","n","phi","v","contamination","innovation","margin","setting","method")
  blocks<-split(seq_len(nrow(d)),interaction(d[groups],drop=TRUE,lex.order=TRUE))
  do.call(rbind,lapply(blocks,function(ix) {
    z<-d[ix,,drop=FALSE]
    out<-z[1L,groups,drop=FALSE]
    out$attempted<-nrow(z);out$finite<-sum(z$valid)
    out$failure_rate<-mean(!z$valid|z$convergence!=0L)
    out$boundary_rate<-if(any(z$valid))mean(z$boundary[z$valid])else NA_real_
    for(p in c("phi","v")) {
      e<-z[[paste0(p,"_hat")]]-z[[p]];e<-e[is.finite(e)];k<-length(e)
      out[[paste0(p,"_bias")]]<-if(k)mean(e)else NA_real_
      out[[paste0(p,"_bias_mcse")]]<-if(k>1L)sd(e)/sqrt(k)else NA_real_
      mse<-if(k)mean(e^2)else NA_real_
      out[[paste0(p,"_rmse")]]<-sqrt(mse)
      out[[paste0(p,"_rmse_mcse")]]<-if(k>1L&&mse>0)sd(e^2)/sqrt(k)/(2*sqrt(mse))else NA_real_
      out[[paste0(p,"_median_abs_error")]]<-if(k)median(abs(e))else NA_real_
    }
    out$spectral_loss<-if(any(is.finite(z$spectral_loss)))mean(z$spectral_loss,na.rm=TRUE)else NA_real_
    out$mean_seconds<-if(any(is.finite(z$elapsed)))mean(z$elapsed,na.rm=TRUE)else NA_real_
    out$mean_fit_seconds<-if(any(is.finite(z$fit_seconds)))mean(z$fit_seconds,na.rm=TRUE)else NA_real_
    out$mean_pilot_seconds<-if(any(is.finite(z$pilot_seconds)))mean(z$pilot_seconds,na.rm=TRUE)else NA_real_
    out$target_interpretation<-z$target_interpretation[1L]
    out$mean_model_calls<-if(any(is.finite(z$calls)))mean(z$calls,na.rm=TRUE)else NA_real_
    out
  }))
}

cse_paired_summary <- function(d) {
  if(is.null(d)||!nrow(d)) return(data.frame())
  key<-c("task","margin","setting")
  base<-d[d$method=="CSE",c(key,"phi_hat","v_hat"),drop=FALSE]
  names(base)[4:5]<-c("phi_cse","v_cse")
  x<-merge(d[d$method!="CSE",,drop=FALSE],base,by=key,all.x=TRUE)
  if(!nrow(x))return(data.frame())
  groups<-c("experiment","n","phi","v","contamination","innovation","margin","setting","method")
  blocks<-split(seq_len(nrow(x)),interaction(x[groups],drop=TRUE,lex.order=TRUE))
  do.call(rbind,lapply(blocks,function(ix) {
    z<-x[ix,,drop=FALSE];o<-z[1L,groups,drop=FALSE]
    for(p in c("phi","v")) {
      a<-z[[paste0(p,"_hat")]]-z[[p]]
      b<-z[[paste0(p,"_cse")]]-z[[p]]
      keep<-is.finite(a)&is.finite(b);a<-a[keep];b<-b[keep]
      delta<-a^2-b^2
      o[[paste0(p,"_paired_n")]]<-length(delta)
      o[[paste0(p,"_mse_comparator_minus_cse")]]<-if(length(delta))mean(delta)else NA_real_
      o[[paste0(p,"_difference_mcse")]]<-if(length(delta)>1L)sd(delta)/sqrt(length(delta))else NA_real_
      o[[paste0(p,"_mse_ratio_to_cse")]]<-if(length(delta)&&mean(b^2)>0)mean(a^2)/mean(b^2)else NA_real_
    }
    o
  }))
}

cse_simulation_variance <- function(d) {
  if(is.null(d)||!any(grepl("CSE_SIM",d$method)))return(data.frame())
  out<-list()
  for(n in unique(d$n))for(phi in unique(d$phi))for(v in unique(d$v)) {
    z<-d[d$n==n&d$phi==phi&d$v==v,,drop=FALSE]
    b<-z[z$method=="CSE"&z$valid,c("phi_hat","v_hat"),drop=FALSE]
    if(nrow(b)<2L)next
    V0<-cov(b)
    for(r in c(1,5,10,20)) {
      a<-z[z$method==paste0("CSE_SIM",r)&z$valid,c("phi_hat","v_hat"),drop=FALSE]
      if(nrow(a)<2L)next
      V<-cov(a)
      out[[length(out)+1L]]<-data.frame(n=n,phi=phi,v=v,ratio=r,
        predicted_factor=1+1/r,phi_variance_ratio=V[1L,1L]/V0[1L,1L],
        v_variance_ratio=V[2L,2L]/V0[2L,2L],complete_cse=nrow(b),complete_sim=nrow(a))
    }
  }
  if(length(out))do.call(rbind,out)else data.frame()
}

cse_write_tex_table <- function(summary,file) {
  if(!nrow(summary))return(invisible(NULL))
  con<-file(file,"w",encoding="UTF-8");on.exit(close(con))
  writeLines(c("\\begingroup", "\\small", "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{longtable}{rrrrllrrrr}","\\toprule",
    "$n$ & $\\phi$ & $v$ & $p$ & Margin & Method & RMSE($\\phi$) & RMSE($v$) & Failure & $N$ \\\\",
    "\\midrule","\\endhead"),con)
  for(i in seq_len(nrow(summary))) {
    z<-summary[i,];method<-gsub("_","\\_",z$method,fixed=TRUE)
    line<-sprintf("%d & %.2f & %.2f & %.3f & %s & %s & %.4f & %.4f & %.3f & %d \\\\",
      z$n,z$phi,z$v,z$contamination,z$margin,method,z$phi_rmse,z$v_rmse,z$failure_rate,z$finite)
    writeLines(line,con)
  }
  writeLines(c("\\bottomrule","\\end{longtable}","\\endgroup"),con)
}

cse_plot_summary <- function(s,file) {
  if(!nrow(s))return(invisible(NULL))
  methods<-unique(s$method);cols<-hcl.colors(length(methods),"Dark 3")
  gr<-unique(s[c("phi","v","contamination","innovation","margin")])
  pdf(file,width=12,height=6.5,onefile=TRUE,family="serif")
  on.exit(dev.off())
  par(mfrow=c(2,4),mar=c(3.8,3.8,2.8,.7),oma=c(0,0,2,0),mgp=c(2.1,.7,0),cex=.85)
  for(i in seq_len(nrow(gr))) {
    g<-gr[i,]
    z<-s[s$phi==g$phi&s$v==g$v&s$contamination==g$contamination&
           s$innovation==g$innovation&s$margin==g$margin,,drop=FALSE]
    for(p in c("phi","v")) {
      yy<-z[[paste0(p,"_rmse")]]
      if(!any(is.finite(yy))) {plot.new();title("No finite estimates");next}
      xr<-range(z$n);if(diff(xr)==0)xr<-xr+c(-.1,.1)*xr[1L]
      ylim<-c(0,max(yy[is.finite(yy)])*1.12+1e-8)
      plot(NA,xlim=xr,ylim=ylim,xlab="Sample size",ylab=paste0("RMSE (",p,")"),
           main=sprintf("phi=%.2f, v=%.2f; %s\np=%.3f, %s",g$phi,g$v,g$margin,g$contamination,g$innovation))
      for(j in seq_along(methods)) {
        a<-z[z$method==methods[j],,drop=FALSE];a<-a[order(a$n),,drop=FALSE]
        if(nrow(a))lines(a$n,a[[paste0(p,"_rmse")]],type="b",pch=j,cex=.7,col=cols[j],lty=j)
      }
      legend("topright",methods,col=cols,lty=seq_along(methods),pch=seq_along(methods),
             cex=.55,bty="n",y.intersp=.8)
    }
  }
}

cse_population <- function(cases,bank,cfg,folder) {
  out<-list();ng<-51L
  axes<-list(phi=seq(cfg$lower[1L],cfg$upper[1L],length.out=ng),v=seq(cfg$lower[2L],cfg$upper[2L],length.out=ng))
  pdf(file.path(folder,"E01_identification.pdf"),width=10,height=8,family="serif")
  on.exit(dev.off());par(mfrow=c(2,2),mar=c(4,4,3,1))
  for(i in seq_len(nrow(cases))) {
    truth<-unlist(cases[i,c("phi","v")],use.names=FALSE);spec<-cse_spec(cases$n[i])
    target<-cse_model(truth,spec,bank)
    Q<-outer(axes$phi,axes$v,Vectorize(function(phi,v)
      sum(spec$weights*(cse_model(c(phi,v),spec,bank)-target)^2)))
    image(axes$phi,axes$v,log10(Q+1e-12),xlab=expression(phi),ylab="v",main="Copula-spectral criterion",col=hcl.colors(32,"YlOrRd",rev=TRUE))
    contour(axes$phi,axes$v,log10(Q+1e-12),add=TRUE,drawlabels=FALSE)
    points(truth[1],truth[2],pch=4,lwd=2)
    image(axes$phi,axes$v,matrix(0,ng,ng),zlim=c(-1,1),xlab=expression(phi),ylab="v",main="Ordinary-spectrum criterion: constant",col="grey90")
    image(axes$phi,axes$v,matrix(0,ng,ng),zlim=c(-1,1),xlab=expression(phi),ylab="v",main="Median-only criterion: constant",col="grey90")
    lambda<-seq(0,pi,length.out=100L)
    ss<-list(tau=c(.1,.2,.5),pairs=rbind(c(1,1),c(1,2),c(3,3)),lambda=lambda)
    M<-matrix(cse_model(truth,ss,bank),100L,3L)
    matplot(lambda/pi,M,type="l",lty=1:3,lwd=1.5,xlab=expression(lambda/pi),ylab="Centered integrated spectrum",main=sprintf("phi=%.2f, v=%.2f",truth[1],truth[2]))
    legend("topright",c("(0.1, 0.1)","(0.1, 0.2)","(0.5, 0.5)"),lty=1:3,col=1:3,bty="n")
    D<-cse_jacobian(truth,spec,bank)
    out[[i]]<-list(truth=truth,axes=axes,criterion=Q,
                   derivative_singular_values=svd(D*sqrt(spec$weights),nu=0,nv=0)$d)
  }
  saveRDS(list(experiment="E01",config=cfg,cases=cases,results=out,numerical_audit=bank$audit),file.path(folder,"E01_results.rds"))
  invisible(out)
}

cse_save_chunk <- function(result,folder,index) {
  path<-file.path(folder,sprintf("chunk_%06d.rds",index))
  if(file.exists(path))stop("Refusing to overwrite a completed checkpoint chunk.")
  temp<-paste0(path,".partial")
  saveRDS(result,temp,compress=FALSE)
  if(!file.rename(temp,path))stop("Could not finalize checkpoint chunk.")
}

cse_run <- function(cfg,source_dir=getwd()) {
  if(!all(cfg$experiments %in% sprintf("E%02d",1:12)))stop("Invalid experiment identifier.")
  if(cfg$cores<1L)stop("CORES must be positive.")
  source_dir<-normalizePath(source_dir,mustWork=TRUE)
  required<-file.path(source_dir,c("cse_core.R","cse_experiments.R"))
  if(!all(file.exists(required)))stop("Both source files must be in source_dir.")
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  cfg$output<-normalizePath(cfg$output,mustWork=TRUE)
  cases<-do.call(rbind,lapply(cfg$experiments,cse_cases,cfg=cfg))
  write.csv(cases,file.path(cfg$output,"experiment_design.csv"),row.names=FALSE)
  tau<-cse_all_tau(cases)
  signature<-list(version=CSE_CODE_VERSION,source_md5=unname(tools::md5sum(required)),tau=tau,
    numerics=cfg[c("lower","upper","trunc_tol","gh_order","bank_nodes","bank_max_nodes","bank_tol","quadrature_tol")])
  cache_dir<-file.path(cfg$output,"checkpoints");dir.create(cache_dir,showWarnings=FALSE)
  cache<-file.path(cache_dir,"model_bank.rds")
  if(file.exists(cache) && cfg$resume) {
    saved<-readRDS(cache)
    if(!identical(saved$signature,signature))stop("Numerical setup changed. Use a new output directory.")
    bank<-saved$bank
  } else {
    message("Building and checking the deterministic model coefficients ...")
    bank_start<-proc.time()[3L]
    bank<-cse_build_bank(tau,cfg)
    bank$audit$setup_seconds<-unname(proc.time()[3L]-bank_start)
    saveRDS(list(signature=signature,bank=bank),cache,compress=FALSE)
  }
  message(sprintf("Model bank ready: J=%d, Gaussian order=%d, variance nodes=%d.",bank$J,bank$audit$gh_order,bank$audit$bank_nodes))
  context<-list(cfg=cfg,bank=bank)
  assign("CSE_CONTEXT",context,envir=.GlobalEnv)
  cl<-NULL
  if(cfg$cores>1L && any(cfg$experiments!="E01")) {
    cl<-parallel::makePSOCKcluster(cfg$cores,outfile="")
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(files) {
      Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
      for(f in files)source(f,local=.GlobalEnv)
      NULL
    },required)
    parallel::clusterExport(cl,"CSE_CONTEXT",envir=.GlobalEnv)
  }
  all_results<-list()
  for(experiment in cfg$experiments) {
    folder<-file.path(cfg$output,experiment);dir.create(folder,showWarnings=FALSE)
    cc<-cases[cases$experiment==experiment,,drop=FALSE]
    if(experiment=="E01") {cse_population(cc,bank,cfg,folder);next}
    reps<-cfg$reps_by_experiment[[experiment]]
    if(is.null(reps))reps<-cfg$reps
    if(length(reps)!=1L || reps<1L)stop("Replication count must be positive.")
    tasks<-cse_task_seeds(cc,cfg,as.integer(reps))
    cp<-file.path(cache_dir,experiment);dir.create(cp,showWarnings=FALSE)
    manifest_path<-file.path(cp,"manifest.rds")
    manifest<-list(signature=signature,cases=cc,reps=reps,
      settings=cfg[setdiff(names(cfg),c("cores","experiments","output","resume","make_plots","checkpoint_every","reps_by_experiment"))])
    if(file.exists(manifest_path) && cfg$resume) {
      if(!identical(readRDS(manifest_path),manifest))stop(paste(experiment,"configuration changed. Use a new output directory."))
    } else {
      if(length(list.files(cp,pattern="^chunk_")))stop("Existing checkpoints found. Use resume=TRUE or a new output directory.")
      saveRDS(manifest,manifest_path)
    }
    chunks<-sort(list.files(cp,pattern="^chunk_[0-9]+\\.rds$",full.names=TRUE))
    results<-if(length(chunks))unlist(lapply(chunks,readRDS),recursive=FALSE)else list()
    done<-vapply(results,function(r)r$id,character(1L))
    pending<-tasks[!vapply(tasks,function(t)t$id %in% done,logical(1L))]
    chunk_index<-length(chunks)
    count<-length(done);total<-length(tasks);start<-proc.time()[3L];new_count<-0L
    message(sprintf("%s: %d/%d tasks already saved; %d workers.",experiment,count,total,cfg$cores))
    if(length(pending)) {
      batches<-split(seq_along(pending),ceiling(seq_along(pending)/max(1L,cfg$checkpoint_every)))
      for(ids in batches) {
        batch<-pending[ids]
        z<-if(is.null(cl))lapply(batch,cse_safe_task)else parallel::parLapplyLB(cl,batch,cse_safe_task)
        chunk_index<-chunk_index+1L;cse_save_chunk(z,cp,chunk_index)
        results<-c(results,z);count<-count+length(z);new_count<-new_count+length(z)
        failed<-sum(vapply(results,function(r)r$status!="completed",logical(1L)))
        elapsed<-unname(proc.time()[3L]-start)
        eta<-if(new_count>0)(total-count)*elapsed/new_count else NA_real_
        message(sprintf("%s %d/%d | task failures %d | elapsed %.1f min | estimated remaining %.1f min",experiment,count,total,failed,elapsed/60,eta/60))
      }
    }
    d<-do.call(rbind,lapply(results,function(r)r$estimates))
    failed_tasks<-lapply(results[vapply(results,function(r)r$status!="completed",logical(1L))],
                         function(r) data.frame(task=r$id,error=r$error,stringsAsFactors=FALSE))
    failure_table<-if(length(failed_tasks))do.call(rbind,failed_tasks)else data.frame(task=character(),error=character())
    s<-cse_summary(d);paired<-cse_paired_summary(d);vr<-cse_simulation_variance(d)
    final<-list(code_version=CSE_CODE_VERSION,config=cfg,cases=cc,
      expected_tasks=length(tasks),completed_tasks=length(results)-nrow(failure_table),
      failed_tasks=failure_table,estimates=d,summary=s,paired=paired,
      simulation_variance=vr,task_results=results,numerical_audit=bank$audit,
      session=utils::sessionInfo(),source_md5=signature$source_md5)
    saveRDS(final,file.path(folder,paste0(experiment,"_results.rds")))
    write.csv(s,file.path(folder,paste0(experiment,"_summary.csv")),row.names=FALSE)
    write.csv(paired,file.path(folder,paste0(experiment,"_paired_comparisons.csv")),row.names=FALSE)
    write.csv(failure_table,file.path(folder,paste0(experiment,"_failed_tasks.csv")),row.names=FALSE)
    if(nrow(vr))write.csv(vr,file.path(folder,paste0(experiment,"_simulation_variance.csv")),row.names=FALSE)
    cse_write_tex_table(s,file.path(folder,paste0(experiment,"_table.tex")))
    if(cfg$make_plots)cse_plot_summary(s,file.path(folder,paste0(experiment,"_RMSE.pdf")))
    all_results[[experiment]]<-list(folder=folder,expected_tasks=length(tasks),failed_tasks=nrow(failure_table))
  }
  saveRDS(list(config=cfg,experiments=all_results,numerical_audit=bank$audit),file.path(cfg$output,"CSE_run_index.rds"))
  message("Finished. Per-experiment RDS files contain estimates, errors, diagnostics, and task RNG states.")
  invisible(all_results)
}
