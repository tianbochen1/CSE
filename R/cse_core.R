CSE_CODE_VERSION <- "1.0.0"

cse_gh <- function(m) {
  A <- matrix(0, m, m)
  A[cbind(seq_len(m - 1L), 2:m)] <- sqrt(seq_len(m - 1L))
  A <- A + t(A)
  e <- eigen(A, symmetric = TRUE)
  o <- order(e$values)
  w <- e$vectors[1L, o]^2
  keep <- w > 0
  list(x = e$values[o][keep], w = w[keep] / sum(w[keep]), order = m)
}

cse_hermite <- function(gh, J) {
  H <- matrix(0, length(gh$x), J + 1L)
  H[, 1L] <- sqrt(gh$w)
  if (J > 0L) H[, 2L] <- gh$x * H[, 1L]
  if (J > 1L) for (j in 1:(J - 1L))
    H[, j + 2L] <- (gh$x * H[, j + 1L] - sqrt(j) * H[, j]) / sqrt(j + 1L)
  H[, -1L, drop = FALSE]
}

cse_quantiles <- function(tau, v, gh) {
  sc <- exp(-sqrt(v) * gh$x / 2 + v / 4)
  vapply(tau, function(u) {
    if (abs(u - 0.5) < 1e-14) return(0)
    p <- min(u, 1 - u)
    fn <- function(q) sum(gh$w * pnorm(q * sc, lower.tail = FALSE)) - p
    hi <- max(1, abs(qnorm(p)))
    while (fn(hi) > 0) hi <- hi * 2
    q <- uniroot(fn, c(0, hi), tol = 1e-11)$root
    if (u < 0.5) -q else q
  }, numeric(1L))
}

cse_coefficients <- function(tau, v, gh, H) {
  q <- cse_quantiles(tau, v, gh)
  sc <- exp(-sqrt(v) * gh$x / 2 + v / 4)
  P <- pnorm(outer(sc, q))
  P <- sweep(P, 2L, tau, "-")
  B <- t(crossprod(H, P * sqrt(gh$w)))
  B[abs(tau - 0.5) < 1e-13, ] <- 0
  list(B = B, q = q)
}

cse_nterms <- function(phi, tol) {
  max(1L, ceiling(log(4 * pi * tol / (1 + 4 * pi * tol)) / log(phi) - 1))
}

cse_interpolate <- function(bank, v, tau, J = bank$J) {
  s <- sqrt(v)
  if (s < bank$s[1L] - 1e-10 || s > tail(bank$s, 1L) + 1e-10)
    stop("Parameter outside the coefficient bank.")
  ix <- match(round(tau, 12L), bank$tau)
  if (anyNA(ix)) stop("Quantile levels are missing from the coefficient bank.")
  i <- min(length(bank$s) - 1L, max(1L, findInterval(s, bank$s)))
  h <- bank$s[i + 1L] - bank$s[i]
  a <- (s - bank$s[i]) / h
  h00 <- 2*a^3 - 3*a^2 + 1
  h10 <- a^3 - 2*a^2 + a
  h01 <- -2*a^3 + 3*a^2
  h11 <- a^3 - a^2
  at <- function(A, k) matrix(A[ix, seq_len(J), k], length(ix), J)
  h00*at(bank$B, i) + h10*h*at(bank$dB, i) +
    h01*at(bank$B, i+1L) + h11*h*at(bank$dB, i+1L)
}

cse_make_pairs <- function(tau, diagonal = FALSE) {
  if (diagonal) return(cbind(seq_along(tau), seq_along(tau)))
  which(upper.tri(matrix(0, length(tau), length(tau)), diag = TRUE), arr.ind = TRUE)
}

cse_spec <- function(n, kind = "full", multiplier = 1L) {
  k <- ceiling(log(n)) * multiplier
  tau <- 0.05 + 0.9 * (seq_len(k) - 0.5) / k
  L <- ceiling(4 * log(n)) * multiplier
  lambda <- pi * (seq_len(L) - 0.5) / L
  if (kind == "central") tau <- 0.3 + 0.4 * (seq_len(k) - 0.5) / k
  if (kind == "tail") {
    a <- 0.05 + 0.15 * (seq_len(ceiling(k/2)) - 0.5) / ceiling(k/2)
    tau <- c(a, 1 - rev(a))
  }
  if (kind == "median") tau <- 0.5
  if (kind == "fixed") {
    tau <- c(0.1, 0.2, 0.35)
    lambda <- pi*c(0.2, 0.5, 0.8)
  }
  if (kind == "low") lambda <- lambda[lambda <= pi/4]
  if (kind == "middle") lambda <- lambda[lambda > pi/4 & lambda < 3*pi/4]
  if (kind == "high") lambda <- lambda[lambda >= 3*pi/4]
  pairs <- cse_make_pairs(tau, kind == "diagonal")
  pw <- ifelse(pairs[, 1L] == pairs[, 2L], 1, 2)
  weights <- rep(pw, each = length(lambda))
  weights <- weights/sum(weights)
  list(tau = round(tau, 12L), lambda = lambda, pairs = pairs,
       weights = weights, kind = kind)
}

cse_model <- function(theta, spec, bank, lags = NULL, direct = FALSE) {
  J <- min(bank$J, cse_nterms(theta[1L], bank$tol))
  B <- if (direct) cse_coefficients(spec$tau, theta[2L], bank$gh,
                                  bank$H[, seq_len(J), drop = FALSE])$B else
    cse_interpolate(bank, theta[2L], spec$tau, J)
  prodB <- B[spec$pairs[,1L],,drop=FALSE] * B[spec$pairs[,2L],,drop=FALSE]
  if (is.null(lags)) {
    r <- theta[1L]^seq_len(J)
    A <- atan(outer(sin(spec$lambda), r) / (1 - outer(cos(spec$lambda), r)))/pi
  } else A <- outer(lags, seq_len(J), function(h,j) theta[1L]^(h*j))
  as.vector(A %*% t(prodB))
}

cse_build_bank <- function(tau, cfg) {
  J <- cse_nterms(cfg$upper[1L], cfg$trunc_tol)
  order <- max(cfg$gh_order, 64L * ceiling((J + 24L)/64L))
  gh <- cse_gh(order)
  H <- cse_hermite(gh, J)
  tau <- sort(unique(round(tau,12L)))
  ns <- cfg$bank_nodes
  repeat {
    s <- seq(sqrt(cfg$lower[2L]), sqrt(cfg$upper[2L]), length.out = ns)
    B <- array(0, c(length(tau), J, ns))
    for (i in seq_len(ns)) B[,,i] <- cse_coefficients(tau,s[i]^2,gh,H)$B
    dB <- B
    step <- s[2L]-s[1L]
    dB[,,1L] <- (-3*B[,,1L]+4*B[,,2L]-B[,,3L])/(2*step)
    dB[,,ns] <- (3*B[,,ns]-4*B[,,ns-1L]+B[,,ns-2L])/(2*step)
    for (i in 2:(ns-1L)) dB[,,i] <- (B[,,i+1L]-B[,,i-1L])/(2*step)
    bank <- list(B=B,dB=dB,s=s,tau=tau,J=J,gh=gh,H=H,tol=cfg$trunc_tol)
    audit_spec <- list(tau=tau, pairs=cse_make_pairs(tau,TRUE),
                       lambda=pi*c(.05,.2,.5,.8,.95))
    vi <- unique(round(seq(1,ns-1L,length.out=min(17L,ns-1L))))
    av <- (s[vi]+0.37*step)^2
    err <- vapply(av,function(v) max(abs(
      cse_model(c(cfg$upper[1L],v),audit_spec,bank)-
      cse_model(c(cfg$upper[1L],v),audit_spec,bank,direct=TRUE))),numeric(1L))
    if (max(err) <= cfg$bank_tol) break
    if (ns >= cfg$bank_max_nodes) stop("Coefficient interpolation audit failed. Increase bank_max_nodes.")
    ns <- min(2L*ns-1L,cfg$bank_max_nodes)
  }
  gh2 <- cse_gh(order+128L)
  H2 <- cse_hermite(gh2,J)
  vv <- c(.03,.25,1,2,2.9)
  vv <- vv[vv>=cfg$lower[2L] & vv<=cfg$upper[2L]]
  qe <- vapply(vv,function(v) {
    b2 <- bank; b2$gh <- gh2; b2$H <- H2
    max(abs(cse_model(c(cfg$upper[1L],v),audit_spec,bank,direct=TRUE)-
              cse_model(c(cfg$upper[1L],v),audit_spec,b2,direct=TRUE)))
  },numeric(1L))
  bank$audit <- list(interpolation=data.frame(v=av,max_error=err),
                     quadrature=data.frame(v=vv,max_error=qe),
                     truncation_bound=cfg$upper[1L]^(J+1L)/(4*pi*(1-cfg$upper[1L]^(J+1L))),
                     gh_order=order,bank_nodes=ns,J=J)
  if (max(qe)>cfg$quadrature_tol)
    stop("Gaussian quadrature audit failed. Increase gh_order and use a new output directory.")
  bank
}

cse_ranks <- function(y) rank(y,ties.method="max")/length(y)

cse_empirical <- function(y, spec, complex = FALSE) {
  n <- length(y); r <- cse_ranks(y)
  A <- outer(r,spec$tau,"<=")*1
  D <- mvfft(A)
  kk <- floor(n*spec$lambda/(2*pi)+1e-12)
  maxk <- max(kk)
  if (!maxk) return(rep(0,length(spec$weights)))
  D <- D[seq_len(maxk)+1L,,drop=FALSE]
  out <- matrix(0+0i,length(kk),nrow(spec$pairs))
  for (p in seq_len(nrow(spec$pairs))) {
    a <- spec$pairs[p,1L]; b <- spec$pairs[p,2L]
    cp <- c(0,cumsum(D[,a]*Conj(D[,b])))/n^2
    g0 <- min(spec$tau[a],spec$tau[b])-spec$tau[a]*spec$tau[b]
    out[,p] <- cp[kk+1L]-spec$lambda*g0/(2*pi)
  }
  if (complex) as.vector(out) else as.vector(Re(out))
}

cse_empirical_lags <- function(y,spec,lags) {
  n <- length(y)
  A <- sweep(outer(cse_ranks(y),spec$tau,"<=")*1,2L,spec$tau,"-")
  ans <- matrix(0,length(lags),nrow(spec$pairs))
  for (i in seq_along(lags)) {
    h <- lags[i]
    M <- crossprod(A[(h+1L):n,,drop=FALSE],A[seq_len(n-h),,drop=FALSE])/(n-h)
    M <- (M+t(M))/2
    ans[i,] <- M[spec$pairs]
  }
  as.vector(ans)
}

cse_innovations <- function(n, distribution="normal") {
  list(initial=rnorm(1L), eta=rnorm(n),
       eps=if(distribution=="t5") rt(n,5)*sqrt(3/5) else rnorm(n))
}

cse_generate <- function(theta,innov) {
  phi <- theta[1L]; v <- theta[2L]; n <- length(innov$eps)
  h <- as.numeric(stats::filter(sqrt(v*(1-phi^2))*innov$eta,
                                filter=phi,method="recursive",init=sqrt(v)*innov$initial))
  exp(h/2-v/4)*innov$eps
}

cse_transform <- function(y,which,v,gh) {
  switch(which, identity=y, asinh=asinh(y), cubic=y+.2*y^3,
         skew=sinh(asinh(y)+.8), cauchy={
           sc <- exp(-sqrt(v)*gh$x/2+v/4)
           surv <- vapply(abs(y),function(z) sum(gh$w*pnorm(z*sc,lower.tail=FALSE)),numeric(1L))
           sign(y)/tan(pi*surv)
         },stop("Unknown transformation."))
}

cse_minimize <- function(fn,cfg,extra=NULL,nonsmooth=FALSE) {
  calls <- 0L
  safe <- function(z) {
    calls <<- calls+1L
    v <- tryCatch(fn(cfg$lower+z*(cfg$upper-cfg$lower)),error=function(e) Inf)
    if(length(v)!=1L || !is.finite(v)) 1e100 else v
  }
  starts <- as.matrix(expand.grid(phi=c(.1,.5,.9),v=c(.1,.5,.9)))
  if (!is.null(extra)) starts <- rbind(starts,
    pmin(1,pmax(0,(extra-cfg$lower)/(cfg$upper-cfg$lower))))
  val <- apply(starts,1L,safe)
  if(all(val>=1e99)) stop("All starting values gave nonfinite criteria.")
  scale <- max(median(val[val<1e99]),1e-9)
  objective <- function(z) safe(z)/scale
  candidates <- order(val)[seq_len(min(cfg$nstarts,length(val)))]
  fits <- vector("list",length(candidates))
  for (i in seq_along(candidates)) {
    z <- starts[candidates[i],]
    if (!nonsmooth) {
      fit <- tryCatch(optim(z,objective,method="L-BFGS-B",lower=c(0,0),upper=c(1,1),
                           control=list(maxit=cfg$maxit,factr=1e7,pgtol=1e-7,ndeps=c(2e-5,2e-5))),
                      error=function(e) NULL)
      if (is.null(fit)) fit <- list(par=z,value=objective(z),convergence=99L)
    } else {
      fit <- cse_pattern(z,objective,cfg$sim_maxit,cfg$sim_parameter_tol)
    }
    fits[[i]] <- fit
  }
  j <- which.min(vapply(fits,function(f) f$value,numeric(1L)))
  best <- fits[[j]]
  z <- pmin(1,pmax(0,best$par))
  list(theta=cfg$lower+z*(cfg$upper-cfg$lower),objective=best$value*scale,
       convergence=best$convergence,calls=calls,boundary=any(z<.002|z>.998))
}

cse_pattern <- function(z,fn,maxit,tol) {
  value <- fn(z); step <- .15; it <- 0L
  directions <- as.matrix(expand.grid(a=-1:1,b=-1:1))
  directions <- directions[rowSums(abs(directions))>0,,drop=FALSE]
  while(step>tol && it<maxit) {
    points <- sweep(step*directions,2L,z,"+")
    points[points<0] <- 0
    points[points>1] <- 1
    vv <- apply(points,1L,fn); j <- which.min(vv)
    if(vv[j] < value-1e-11*(1+abs(value))) {z<-points[j,];value<-vv[j]}
    else step <- step/2
    it <- it+1L
  }
  list(par=z,value=value,convergence=if(step<=tol) 0L else 1L)
}

cse_fit_spectral <- function(y,spec,bank,cfg,W=NULL,extra=NULL) {
  empirical <- cse_empirical(y,spec)
  fn <- function(theta) {
    d <- empirical-cse_model(theta,spec,bank)
    if(is.null(W)) sum(spec$weights*d*d) else drop(crossprod(d,W%*%d))
  }
  cse_minimize(fn,cfg,extra)
}

cse_fit_td <- function(y,spec,lags,bank,cfg) {
  empirical <- cse_empirical_lags(y,spec,lags)
  pw <- ifelse(spec$pairs[,1L]==spec$pairs[,2L],1,2)
  weights <- rep(pw,each=length(lags));weights<-weights/sum(weights)
  fn <- function(theta) sum(weights*(empirical-cse_model(theta,spec,bank,lags))^2)
  cse_minimize(fn,cfg)
}

cse_qml_value <- function(theta,x) {
  phi<-theta[1L];v<-theta[2L];R<-pi^2/2;Q<-v*(1-phi^2)
  a<-0;b<-0;P<-v;srr<-0;srd<-0;sdd<-0;ld<-0
  for (xt in x) {
    F<-P+R;r<-xt-a;d<-1+b
    srr<-srr+r*r/F;srd<-srd+r*d/F;sdd<-sdd+d*d/F;ld<-ld+log(F)
    K<-P/F;a<-phi*(a+K*r);b<-phi*(b-K*d);P<-phi^2*(P-K*P)+Q
  }
  (ld+srr-srd*srd/sdd)/length(x)
}

cse_fit_qml <- function(y,cfg) {
  if(any(!is.finite(y)) || any(y==0)) stop("Log-squared QML requires finite nonzero observations.")
  x<-2*log(abs(y))
  cse_minimize(function(theta)cse_qml_value(theta,x),cfg)
}

cse_fit_simulated <- function(y,spec,cfg,innov,extra) {
  empirical <- cse_empirical(y,spec,complex=TRUE)
  fn <- function(theta) {
    model <- cse_empirical(cse_generate(theta,innov),spec,complex=TRUE)
    sum(spec$weights*Mod(empirical-model)^2)
  }
  cse_minimize(fn,cfg,extra,nonsmooth=TRUE)
}

cse_covariance_weight <- function(theta,n,spec,B,cfg) {
  p<-length(spec$weights)
  draws<-matrix(0,B,p)
  for(i in seq_len(B)) draws[i,]<-cse_empirical(cse_generate(theta,cse_innovations(n)),spec)
  S<-cov(draws)
  ee<-eigen(S,symmetric=TRUE)
  floor<-max(max(ee$values)*cfg$weight_ridge,1e-12)
  W<-ee$vectors%*%(t(ee$vectors)/pmax(ee$values,floor))
  W<-W/mean(diag(W))
  list(W=W,eigenvalues=ee$values,ridge=floor)
}

cse_cpl_grid <- function(n,m) plogis(seq(qlogis(1/(n+1)),qlogis(n/(n+1)),length.out=m))

cse_bilinear_setup <- function(u1,u2,grid) {
  i<-pmin(length(grid)-1L,pmax(1L,findInterval(u1,grid)))
  j<-pmin(length(grid)-1L,pmax(1L,findInterval(u2,grid)))
  list(i=i,j=j,a=(u1-grid[i])/(grid[i+1L]-grid[i]),b=(u2-grid[j])/(grid[j+1L]-grid[j]))
}

cse_bilinear <- function(M,s) {
  (1-s$a)*(1-s$b)*M[cbind(s$i,s$j)]+s$a*(1-s$b)*M[cbind(s$i+1L,s$j)]+
    (1-s$a)*s$b*M[cbind(s$i,s$j+1L)]+s$a*s$b*M[cbind(s$i+1L,s$j+1L)]
}

cse_logcopula_grid <- function(theta,grid,lags,gh,gh_margin) {
  v<-theta[2L];q<-cse_quantiles(grid,v,gh_margin)
  scale<-exp(-sqrt(v)*gh$x/2+v/4)
  A<-dnorm(outer(q,scale))*rep(scale,each=length(q))
  ms<-exp(-sqrt(v)*gh_margin$x/2+v/4)
  marginal<-drop((dnorm(outer(q,ms))*rep(ms,each=length(q)))%*%gh_margin$w)
  if(any(marginal<=0)) stop("Nonpositive marginal density in copula likelihood.")
  lapply(lags,function(h) {
    r<-theta[1L]^h
    C<-matrix(0,length(q),length(gh$x))
    for(i in seq_along(gh$x)) {
      cond<-r*gh$x[i]+sqrt(1-r*r)*gh$x
      sc<-exp(-sqrt(v)*cond/2+v/4)
      C[,i]<-drop((dnorm(outer(q,sc))*rep(sc,each=length(q)))%*%gh$w)
    }
    joint<-(A*rep(gh$w,each=length(q)))%*%t(C)
    joint<-(joint+t(joint))/2
    if(any(joint<=0) || any(!is.finite(joint))) stop("Invalid copula density quadrature.")
    log(joint)-outer(log(marginal),log(marginal),"+")
  })
}

cse_fit_cpl <- function(y,bank,cfg,extra=NULL,refined=FALSE) {
  n<-length(y);u<-rank(y,ties.method="max")/(n+1)
  m<-if(refined) 2L*cfg$cpl_grid-1L else cfg$cpl_grid
  grid<-cse_cpl_grid(n,m)
  gh<-cse_gh(if(refined) cfg$cpl_gh+16L else cfg$cpl_gh)
  lags<-c(1L,2L)
  interp<-lapply(lags,function(h)cse_bilinear_setup(u[(h+1L):n],u[seq_len(n-h)],grid))
  fn<-function(theta) {
    M<-cse_logcopula_grid(theta,grid,lags,gh,bank$gh)
    -sum(vapply(seq_along(lags),function(k)sum(cse_bilinear(M[[k]],interp[[k]])),numeric(1L)))/sum(n-lags)
  }
  ans<-cse_minimize(fn,cfg,extra)
  ans
}

cse_jacobian <- function(theta,spec,bank) {
  D<-matrix(0,length(spec$weights),2L)
  for(j in 1:2) {
    delta<-1e-5; a<-b<-theta;a[j]<-a[j]+delta;b[j]<-b[j]-delta
    D[,j]<-(cse_model(a,spec,bank)-cse_model(b,spec,bank))/(2*delta)
  }
  D
}
