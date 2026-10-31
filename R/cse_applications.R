cse_sources <- function() {
  data.frame(application=c("EURUSD","SP500"),series=c("DEXUSEU","SP500"),
    label=c("EUR/USD exchange rate","S&P 500 price index"),
    provider=c("Board of Governors of the Federal Reserve System","S&P Dow Jones Indices LLC"),
    unit=c("US dollars per euro","Price index; dividends excluded"),stringsAsFactors=FALSE)
}

cse_read_prices <- function(file, series) {
  d <- read.csv(file,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c(".","NA",""))
  nm <- tolower(names(d)); di <- which(nm %in% c("date","observation_date"))[1L]
  vi <- match(series,names(d))
  if (is.na(vi)) vi <- which(nm %in% c("value","price","close"))[1L]
  if (is.na(di) || is.na(vi)) stop("CSV needs a date/observation_date column and ",series,"/value column.")
  ans <- data.frame(date=as.Date(d[[di]]),value=suppressWarnings(as.numeric(d[[vi]])))
  if (anyNA(ans$date)) stop("Invalid dates in ",file)
  if (anyDuplicated(ans$date)) stop("Duplicate dates in ",file)
  ans <- ans[order(ans$date),]
  if (!nrow(ans) || sum(is.finite(ans$value))<100) stop("Too few valid prices in ",file)
  ans
}

cse_load_application <- function(app,cfg) {
  meta <- cse_sources(); meta <- meta[meta$application==app,,drop=FALSE]
  dir.create(cfg$data_dir,recursive=TRUE,showWarnings=FALSE)
  file <- file.path(cfg$data_dir,paste0(meta$series,".csv"))
  local <- cfg$local_csv[[app]]
  if (!is.null(local) && nzchar(local)) {
    if (!file.exists(local)) stop("LOCAL_CSV file does not exist: ",local)
    file <- normalizePath(local,mustWork=TRUE)
  } else if (cfg$refresh || !file.exists(file)) {
    u <- paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=",meta$series,
                "&cosd=",cfg$start-21L,"&coed=",cfg$end)
    tmp <- tempfile(fileext=".csv",tmpdir=cfg$data_dir)
    on.exit(unlink(tmp),add=TRUE)
    old <- options(timeout=max(60,getOption("timeout",60))); on.exit(options(old),add=TRUE)
    error <- NULL
    for (attempt in 1:2) {
      ok <- tryCatch({status<-download.file(u,tmp,mode="wb",quiet=TRUE,method="libcurl")
        if(status!=0L) stop("download.file returned ",status)
        cse_read_prices(tmp,meta$series); TRUE},error=function(e){error<<-conditionMessage(e);FALSE})
      if(ok) break
    }
    if(!ok) stop("Download failed for ",app,": ",error,
      ". Download the FRED CSV manually and set LOCAL_CSV$",app," to its path.")
    if(!file.copy(tmp,file,overwrite=TRUE)) stop("Cannot save data cache.")
  }
  raw <- cse_read_prices(file,meta$series)
  requested <- raw$date>=cfg$start & raw$date<=cfg$end
  missing_prices <- sum(!is.finite(raw$value[requested]))
  prices <- raw[is.finite(raw$value),]
  prices <- prices[prices$date>=cfg$start-21L & prices$date<=cfg$end,]
  if (min(prices$date)>cfg$start+7L || max(prices$date)<cfg$end-7L)
    stop("Cached/downloaded data do not cover the requested period. Refresh the cache or supply a complete local CSV.")
  if(any(prices$value<=0)) stop("Nonpositive price: log returns cannot be computed for ",app)
  d <- data.frame(date=prices$date[-1L],previous_date=prices$date[-nrow(prices)],
                  price=prices$value[-1L],y=100*diff(log(prices$value)))
  d <- d[d$date>=cfg$start & d$date<=cfg$end,]; rownames(d)<-NULL
  if(nrow(d)<max(512L,cfg$train_n+cfg$test_n)) stop("Insufficient observations for full fit and validation: ",app)
  gap <- as.integer(d$date-d$previous_date)
  if(max(gap)>10L) stop("An observation gap exceeds 10 calendar days; inspect missing data before analysis.")
  freq <- table(d$y)
  audit <- data.frame(application=app,series=meta$series,start=min(d$date),end=max(d$date),
       n=nrow(d),missing_price_rows=missing_prices,zero_returns=sum(d$y==0),
       duplicate_return_fraction=1-length(unique(d$y))/nrow(d),
       largest_tie_fraction=max(freq)/nrow(d),max_calendar_gap=max(gap),
       mean_return=mean(d$y),sd_return=sd(d$y),
       source_md5=unname(tools::md5sum(file)),stringsAsFactors=FALSE)
  list(meta=meta,data=d,audit=audit,source_file=file,
       source_url=paste0("https://fred.stlouisfed.org/series/",meta$series),
       retrieved_on=as.character(Sys.Date()))
}

cse_eval_spec <- function() {
  tau <- c(.1,.2,.35,.5,.65,.8,.9); lambda<-pi*(seq_len(32L)-.5)/32L
  pairs<-cse_make_pairs(tau); w<-rep(ifelse(pairs[,1]==pairs[,2],1,2),each=length(lambda))
  list(tau=tau,lambda=lambda,pairs=pairs,weights=w/sum(w),kind="validation")
}

cse_fit_real <- function(y,method,spec,bank,cfg,extra=NULL) {
  started<-proc.time()[[3L]]
  ans<-switch(method,
    CSE=cse_fit_spectral(y,spec,bank,cfg,extra=extra),
    TD2=cse_fit_td(y,spec,c(1L,2L),bank,cfg),
    TD16=cse_fit_td(y,spec,c(1L,2L,4L,8L,16L),bank,cfg),
    QML=cse_fit_qml(y-mean(y),cfg),
    CPL=cse_fit_cpl(y,bank,cfg,extra=extra,refined=TRUE),
    stop("Unknown method: ",method))
  ans$elapsed<-unname(proc.time()[[3L]]-started); ans$method<-method
  if(any(!is.finite(ans$theta))) stop("Nonfinite estimates for ",method)
  ans
}

cse_validation_score <- function(y,theta,spec,bank) {
  e<-cse_empirical(y,spec,complex=TRUE)
  model<-if(is.null(theta)) rep(0,length(e)) else cse_model(theta,spec,bank)
  lags<-seq_len(20L); el<-cse_empirical_lags(y,spec,lags)
  ml<-if(is.null(theta)) rep(0,length(el)) else cse_model(theta,spec,bank,lags)
  w<-rep(ifelse(spec$pairs[,1]==spec$pairs[,2],1,2),each=length(lags));w<-w/sum(w)
  c(spectral_real=sum(spec$weights*(Re(e)-model)^2),
    spectral_complex=sum(spec$weights*Mod(e-model)^2),
    lag_discrepancy=sum(w*(el-ml)^2))
}

cse_check_stats <- function(y,theta,spec,bank) {
  empirical<-cse_empirical(y,spec,complex=TRUE); m<-cse_model(theta,spec,bank)
  a<-as.numeric(cse_ranks(y)<=.5)
  ac<-as.numeric(acf(a,lag.max=20L,plot=FALSE)$acf)[-1L]
  list(statistics=c(global=length(y)*sum(spec$weights*Mod(empirical-m)^2),
        imaginary=length(y)*sum(spec$weights*Im(empirical)^2),median=length(y)*sum(ac^2)),
       empirical=empirical)
}

cse_rows <- function(xs) {
  if(!length(xs)) return(data.frame())
  out<-do.call(rbind,xs);rownames(out)<-NULL;out
}

cse_fit_row <- function(fit,method=fit$method) {
  data.frame(method=method,phi=fit$theta[1L],v=fit$theta[2L],
      latent_half_life=log(.5)/log(fit$theta[1L]),objective=fit$objective,
      boundary=fit$boundary,optimizer_code=fit$convergence,seconds=fit$elapsed,
      stringsAsFactors=FALSE)
}


cse_application_spec <- function(n, q_mult=1L, f_mult=1L) {
  k<-ceiling(log(n))*q_mult; L<-ceiling(4*log(n))*f_mult
  tau<-round(.05+.9*(seq_len(k)-.5)/k,12L)
  lambda<-pi*(seq_len(L)-.5)/L; pairs<-cse_make_pairs(tau)
  w<-rep(ifelse(pairs[,1L]==pairs[,2L],1,2),each=L)
  list(tau=tau,lambda=lambda,pairs=pairs,weights=w/sum(w),kind="full")
}


cse_component_checks <- function(y,theta,spec,bank) {
  z<-cse_check_stats(y,theta,spec,bank)
  z$statistics<-c(global=z$statistics[["global"]],
    real=z$statistics[["global"]]-z$statistics[["imaginary"]],
    imaginary=z$statistics[["imaginary"]],median=z$statistics[["median"]])
  z
}


cse_model <- function(theta,spec,bank,lags=NULL,direct=FALSE) {
  J <- min(bank$J,cse_nterms(theta[1L],bank$tol))
  B <- if(direct) cse_coefficients(spec$tau,theta[2L],bank$gh,bank$H[,seq_len(J),drop=FALSE])$B else cse_interpolate(bank,theta[2L],spec$tau,J)
  # Bound the contribution of discarded Hermite terms by 1e-13.
  r <- theta[1L]^seq_len(J)
  envelopes <- apply(abs(B),2L,max)^2
  envelopes <- envelopes * if(is.null(lags)) r/(pi*(1-r)) else r^min(lags)
  remaining <- c(rev(cumsum(rev(envelopes)))[-1L],0)
  keep <- which(remaining<=1e-13)[1L]
  B <- B[,seq_len(keep),drop=FALSE];r<-r[seq_len(keep)]
  prodB <- B[spec$pairs[,1L],,drop=FALSE]*B[spec$pairs[,2L],,drop=FALSE]
  if(is.null(lags)) A<-atan(outer(sin(spec$lambda),r)/(1-outer(cos(spec$lambda),r)))/pi else
    A<-outer(lags,seq_len(keep),function(h,j)theta[1L]^(h*j))
  as.vector(A%*%t(prodB))
}

cse_objective_spectral <- function(y,spec,bank) {
  n<-length(y);k<-length(spec$tau);L<-length(spec$lambda)
  stopifnot(nrow(spec$pairs)==k*(k+1)/2)
  empirical<-cse_empirical(y,spec)
  constant<-sum(spec$weights*empirical^2)
  rm(empirical)
  kk<-floor(n*spec$lambda/(2*pi)+1e-12)
  d<-mvfft(outer(cse_ranks(y),spec$tau,"<=")*1)[seq_len(max(kk))+1L,,drop=FALSE]
  dr<-Re(d);di<-Im(d)
  g0<-outer(spec$tau,spec$tau,pmin)-outer(spec$tau,spec$tau)
  spec$weights<-NULL;spec$pairs<-NULL
  gc(FALSE)
  function(theta) {
    J<-min(bank$J,cse_nterms(theta[1],bank$tol))
    B<-cse_interpolate(bank,theta[2],spec$tau,J)
    r<-theta[1]^seq_len(J)
    e<-apply(abs(B),2,max)^2*r/(pi*(1-r))
    tail<-c(rev(cumsum(rev(e)))[-1],0)
    keep<-which(tail<=1e-13)[1]
    B<-B[,seq_len(keep),drop=FALSE];r<-r[seq_len(keep)]
    angle<-atan(outer(sin(spec$lambda),r)/(1-outer(cos(spec$lambda),r)))/pi
    z<-(dr%*%B)^2+(di%*%B)^2
    z<-rbind(rep(0,keep),apply(z,2,cumsum))
    h<-z[kk+1L,,drop=FALSE]/n^2-outer(spec$lambda/(2*pi),colSums(B*(g0%*%B)))
    gram<-crossprod(B)^2
    constant+(sum((angle%*%gram)*angle)-2*sum(angle*h))/(L*k^2)
  }
}
cse_fit_spectral <- function(y,spec,bank,cfg,W=NULL,extra=NULL) {
  if(!is.null(W))stop("This accelerated criterion requires the full uniformly weighted grid.")
  cse_minimize(cse_objective_spectral(y,spec,bank),cfg,extra)
}

