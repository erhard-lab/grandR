logSumExp <- function(log_values) {
  max_log <- max(log_values)
  max_log + log(sum(exp(log_values - max_log)))
}
fit.ntr.betamix <- function(ll0,ll1,c,
                            beta.approx = TRUE,
                            conversion_reads = FALSE,
                            len = 100,nstart = 20,plot=FALSE) {
  start <- 0
  optfun = function(p) {
    ldb = lse(ll0+log(1-p),ll1+log(p))
    sum(c * ldb) - start
  }
  opt <- optimize(optfun, maximum = TRUE, lower = 0, upper = 1)

  ll_0 <- optfun(0)
  ll_1 <- optfun(1)

  if (ll_0 > opt$objective) {
    opt$maximum <- 0
    opt$objective <- ll_0
  } else if (ll_1 > opt$objective) {
    opt$maximum <- 1
    opt$objective <- ll_1
  }
  start <- opt$objective + start

  left <- if (optfun(0) > log(1E-3)) 0 else uniroot(function(x) optfun(x) - log(1E-3), c(0, opt$maximum))$root
  right <- if (optfun(1) > log(1E-3)) 1 else uniroot(function(x) optfun(x) - log(1E-3), c(opt$maximum, 1))$root

  if (!beta.approx) {
    return(if (conversion.reads)
      c(ntr = opt$maximum, conversion.reads = sum(mixmat) - sum(mixmat[0, ]))
      else opt$maximum)
  }
  if ((opt$maximum == 0  && (left == 0 && right == 1))) { x <- qbeta(seq(left, right, length.out = len), 0.5, 0.9)
  } else if ((opt$maximum == 1 && (left == 0 && right == 1))) { x <- qbeta(seq(left, right, length.out = len), 0.9, 0.5)
  } else {
    x = seq(left, right, length.out = len)
  }
  logpost <- sapply(x, optfun)
  log_int <- logSumExp(logpost) + log(diff(range(x)) / (len - 1))
  integral <- start + log_int
  log_density <- logpost - logSumExp(logpost)
  density <- exp(log_density)
  dx <- diff(x)
  mid_vals <- (density[-1] + density[-length(density)]) / 2
  fs2 <- c(0, cumsum(mid_vals * dx))
  fs2 <- fs2 / max(fs2)
  fit <- fit_beta_mixture_cdf_fast(x, fs2,opt$maximum, nstart = nstart)

  # dx <- diff(x)
  # fx_obs <- (fs2[-1] + fs2[-length(fs2)]) / 2
  # fx_fit <- (mix_cdf(x[-1]) + mix_cdf(x[-length(x)])) / 2
  # error_mix <- sum((fx_obs - fx_fit)^2 * dx)
  if (plot) {
    mix_cdf <- function(xnew) {
      fit["w"] * pbeta(xnew, fit["a1"], fit["b1"]) +
        (1 - fit["w"]) * pbeta(xnew, fit["a2"], fit["b2"])
    }
    plot(x,fs2,xlab="ntr",ylab="Cumulative freq",type='l')
    graphics::lines(x,mix_cdf(x),col='red')
    graphics::legend("topleft",legend=c("Actual distribution","Beta Mix approximation"),fill=c("black","red"))
  }
  return(c(ntr = opt$maximum, fit, int = integral))
}

fit_beta_mixture_cdf_fast <- function(x, Femp,
                                      opt_maximum = NULL,
                                      nstart      = 5,
                                      p_range     = c(0,1),
                                      shape_lower = 1) {

  stopifnot(length(x) == length(Femp), all(diff(x) > 0))

  dx   <- diff(x)
  xm   <- (head(x, -1) + tail(x, -1)) / 2
  dens <- pmax(diff(Femp) / dx, 1e-12)
  weights <- dens * dx
  weights <- weights / sum(weights)

  sqErr <- function(par) {
    Fmix <- par[1] * pbeta(x, par[2], par[3]) +
      (1 - par[1]) * pbeta(x, par[4], par[5])
    sum((Femp - Fmix)^2)
  }

  lower <- c(p_range[1], rep(shape_lower, 4))
  upper <- c(p_range[2], rep(Inf, 4))

  run_optim <- function(init) {
    tryCatch(
      optim(init, sqErr, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 1000)),
      error = function(e) NULL
    )
  }

  single_sqErr <- function(par) {
    Fbeta <- pbeta(x, par[1], par[2])
    sum((Femp - Fbeta)^2)
  }

  sf <- optim(c(3, 3), single_sqErr, method = "L-BFGS-B", lower = c(shape_lower, shape_lower))
  single_fit = sf$par
  Fsingle <- pbeta(x, single_fit[1], single_fit[2])
  single_SSE <- sum((Femp - Fsingle)^2)

  if (single_SSE > 1e-7) {

    p_init <- 0.7
    ab_init <- rbind(c(pmax(single_fit[1],shape_lower+0.1), pmax(single_fit[2],shape_lower+0.1)), c(2, 2))
    mom_init <- c(p_init, as.vector(t(ab_init)))
    res <- run_optim(mom_init)

    if (is.null(res) || res$convergence != 0 || res$value >= single_SSE * (1 - 0.5)) {
      multi_start_attempted <- TRUE

      for (attempt in 1:nstart) {
        random_init <- c(runif(1, p_range[1], p_range[2]),
                         runif(4, shape_lower, 5))
        res <- run_optim(random_init)
        attempts_used <- attempt + 1
        if (!is.null(res) && res$convergence == 0 && res$value <= single_SSE * (1 - 0.5)) break
      }
      m1 <- res$par[2] / (res$par[2] + res$par[3])
      m2 <- res$par[4] / (res$par[4] + res$par[5])

      if (m1 > m2) {
        res$par <- c(1 - res$par[1], res$par[4], res$par[5], res$par[2], res$par[3])
      }
    }

    if (!is.null(res) && res$convergence == 0 && res$value <= single_SSE) {
      return(setNames(c(res$par,
                        SSE = res$value),
                      c("w", "a1", "b1", "a2", "b2", "SSE")))

    }
  }

  return(setNames(c(w = 1, a1 = single_fit[1], b1 = single_fit[2], a2 = 2, b2 = 2,
                    SSE = single_SSE),
                  c("w", "a1", "b1", "a2", "b2", "SSE")))
}
