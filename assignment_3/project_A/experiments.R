# Source from the project root, or source files explicitly before using these functions.
fit_all <- function(x, nu, start=c(median(x),log(IQR(x)^2/4)), tol=1e-7) {
  list(EM=fit_em(x,nu,start,tol), BFGS=fit_bfgs(x,nu,start,tol), Newton=fit_newton(x,nu,start,tol))
}
fit_table <- function(fits) {
  do.call(rbind,lapply(names(fits),function(method) {
    f <- fits[[method]]
    data.frame(method,mu=f$mu,v=f$v,loglik=f$loglik,iterations=f$iterations,
               objective=f$evaluations[1],gradient=f$evaluations[2],hessian=f$evaluations[3],
               converged=f$converged,seconds=f$elapsed,row.names=NULL)
  }))
}
# Independent reference: derivative-free Nelder-Mead from several starting points.
reference_optimum <- function(x, nu) {
  starts <- list(c(median(x),log(IQR(x)^2/4)),c(mean(x),log(var(x))),c(median(x)+1,0))
  values <- sapply(starts,function(s) optim(s,function(e) -t_loglik(x,e[1],exp(e[2]),nu),
                          method="Nelder-Mead",control=list(reltol=1e-13,maxit=20000))$value)
  -min(values)
}
# Repeated timing; accuracy-certified fits only are eligible for runtime ranking.
benchmark_methods <- function(n=1000,repetitions=7,tol=1e-7) {
  rows <- list(); index <- 0
  for(nu in c(2,5,30)) {
    x <- simulate_t_data(n,1,2,nu,seed=400+nu)
    fit_all(x,nu,tol=tol) # warm up before collecting timings
    for(repetition in seq_len(repetitions)) {
      fits <- fit_all(x,nu,tol=tol)
      index <- index+1; rows[[index]] <- cbind(nu,repetition,fit_table(fits))
    }
  }
  do.call(rbind,rows)
}
# Sample the readable implementation, then compare a targeted validation improvement.
profile_em <- function(n=20000,repetitions=25,interval=0.002) {
  x <- simulate_t_data(n,1,2,2,seed=188)
  path <- tempfile(); on.exit(unlink(path))
  Rprof(path,interval=interval)
  tryCatch(for(i in seq_len(repetitions)) fit_em(x,2),finally=Rprof(NULL))
  profile <- summaryRprof(path)
  rows <- list()
  for(i in seq_len(7)) for(fast in c(FALSE,TRUE)) {
    f <- fit_em(x,2,optimized=fast)
    rows[[length(rows)+1]] <- data.frame(optimized=fast,seconds=f$elapsed,mu=f$mu,v=f$v,loglik=f$loglik,converged=f$converged)
  }
  list(profile=profile,comparison=do.call(rbind,rows),n=n,repetitions=repetitions)
}

# Monte Carlo study: paired datasets, two starts, and all three optimizers.
# Returns individual fit records and per-dataset information-based uncertainty.
# An optional RDS cache is invalidated when settings or source files change.
run_monte_carlo <- function(repetitions=300, n=300, nu_values=c(2,5,30),
                            seed=32026, tol=1e-7, maxit=200,
                            cache_file=NULL, source_dir="project_A", force=FALSE) {
  settings <- list(repetitions=repetitions,n=n,nu_values=nu_values,seed=seed,tol=tol,maxit=maxit)
  if (repetitions < 2 || n < 3) stop("Use at least two replications and three observations")
  source_files <- sort(list.files(source_dir,pattern="\\.R$",full.names=TRUE))
  signature <- list(settings=settings,source_md5=unname(tools::md5sum(source_files)),R=R.version.string)
  if (!is.null(cache_file) && file.exists(cache_file) && !force) {
    saved <- readRDS(cache_file)
    if (identical(saved$signature,signature)) return(saved)
  }
  methods <- list(EM=fit_em,BFGS=fit_bfgs,Newton=fit_newton)
  fit_rows <- list(); uncertainty_rows <- list(); agreement_rows <- list()
  # Warm-up uses a separate seed and is excluded from recorded fit times.
  warmup <- simulate_t_data(n,0,1,5,seed=seed)
  invisible(fit_all(warmup,5,tol=tol))
  for (j in seq_along(nu_values)) {
    nu <- nu_values[j]
    for (replication in seq_len(repetitions)) {
      dataset_seed <- seed+1000*j+replication
      x <- simulate_t_data(n,0,1,nu,seed=dataset_seed)
      initial_scale <- IQR(x)/2
      starts <- list("Data-based"=c(median(x),2*log(initial_scale)),
                     "Displaced"=c(median(x)+5*initial_scale,2*log(initial_scale/2)))
      for (scenario in names(starts)) {
        fits <- list()
        # Rotate execution order to reduce systematic first/last-method timing effects.
        order <- ((seq_along(methods)+replication-2) %% length(methods))+1
        for (method in names(methods)[order]) {
          tic <- proc.time()[3]
          fit <- tryCatch(methods[[method]](x,nu,start=starts[[scenario]],tol=tol,maxit=maxit),
                          error=function(e) e)
          seconds <- unname(proc.time()[3]-tic)
          if (inherits(fit,"error")) {
            row <- data.frame(mu=NA_real_,v=NA_real_,loglik=NA_real_,converged=FALSE,
                              score=NA_real_,reason=conditionMessage(fit))
          } else {
            error <- score_error(x,c(fit$mu,log(fit$v)),nu)
            success <- isTRUE(fit$converged) && is.finite(error) && error<=tol
            row <- data.frame(mu=fit$mu,v=fit$v,loglik=fit$loglik,converged=success,score=error,
                              reason=if(success) "Certified" else "Stopping rule / score not certified")
            if (success) fits[[method]] <- fit
          }
          fit_rows[[length(fit_rows)+1]] <- cbind(data.frame(nu,replication,dataset_seed,scenario,method,seconds),row)
        }
        # Sampling uncertainty is studied once per dataset, using the data-based start.
        # Choose the highest likelihood among certified fits; never count methods as replicates.
        if (scenario=="Data-based") {
          agreement_rows[[length(agreement_rows)+1]] <- data.frame(nu,replication,
            certified=length(fits),
            mu_spread=if(length(fits)>1) diff(range(sapply(fits,`[[`,"mu"))) else NA_real_,
            sigma_spread=if(length(fits)>1) diff(range(sqrt(sapply(fits,`[[`,"v")))) else NA_real_,
            loglik_spread=if(length(fits)>1) diff(range(sapply(fits,`[[`,"loglik"))) else NA_real_)
          if (length(fits)==0) next
          fit <- fits[[which.max(sapply(fits,`[[`,"loglik"))]]
          for (information_type in c("Expected","Observed")) {
            information <- if(information_type=="Expected") expected_information(fit$v,nu,n) else
              observed_information(x,fit$mu,fit$v,nu)
            covariance <- tryCatch({chol(information); solve(information)},error=function(e) NULL)
            # Delta method from (mu,v) to (mu,sigma): d sigma / d v = 1/(2 sqrt(v)).
            standard_errors <- if(is.null(covariance)) c(NA_real_,NA_real_) else
              sqrt(diag(covariance))*c(1,1/(2*sqrt(fit$v)))
            estimates <- c(fit$mu,sqrt(fit$v)); truth <- c(0,1)
            for (parameter in 1:2) {
              se <- standard_errors[parameter]
              uncertainty_rows[[length(uncertainty_rows)+1]] <- data.frame(nu,replication,
                information=information_type,parameter=c("mu","sigma")[parameter],
                estimate=estimates[parameter],truth=truth[parameter],se=se,
                covered=if(is.finite(se)) abs(estimates[parameter]-truth[parameter])<=qnorm(.975)*se else NA)
            }
          }
        }
      }
    }
  }
  result <- list(signature=signature,settings=settings,fits=do.call(rbind,fit_rows),
                 uncertainty=do.call(rbind,uncertainty_rows),agreement=do.call(rbind,agreement_rows),
                 completed_at=format(Sys.time(),tz="Europe/Copenhagen",usetz=TRUE),session=sessionInfo())
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file),recursive=TRUE,showWarnings=FALSE)
    saveRDS(result,cache_file)
  }
  result
}

# Compact calibration table. Missing information inversions are counted explicitly.
uncertainty_summary <- function(study) {
  groups <- split(study$uncertainty,interaction(study$uncertainty[c("nu","parameter","information")],drop=TRUE))
  rows <- lapply(groups,function(group) {
    usable <- is.finite(group$se)
    empirical_sd <- sd(group$estimate)
    average_se <- if(any(usable)) mean(group$se[usable]) else NA_real_
    coverage <- if(any(usable)) mean(group$covered[usable]) else NA_real_
    data.frame(nu=group$nu[1],parameter=group$parameter[1],information=group$information[1],
      datasets=nrow(group),valid_SE=sum(usable),bias=mean(group$estimate-group$truth),
      empirical_SD=empirical_sd,average_SE=average_se,SE_SD_ratio=average_se/empirical_sd,
      SD_MCSE=empirical_sd/sqrt(2*(nrow(group)-1)),coverage=coverage,
      coverage_MCSE=if(any(usable)) sqrt(coverage*(1-coverage)/sum(usable)) else NA_real_)
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result[order(result$parameter,result$nu,result$information),]
}

# Wilson intervals quantify uncertainty in failure rates, including zero failures.
# Timing summaries include failures: failed fits still consume computational time.
runtime_summary <- function(study) {
  groups <- split(study$fits,interaction(study$fits[c("nu","scenario","method")],drop=TRUE))
  rows <- lapply(groups,function(group) {
    total <- nrow(group); failures <- sum(!group$converged); p <- failures/total; z <- qnorm(.975)
    center <- (p+z^2/(2*total))/(1+z^2/total)
    half <- z*sqrt(p*(1-p)/total+z^2/(4*total^2))/(1+z^2/total)
    data.frame(nu=group$nu[1],scenario=group$scenario[1],method=group$method[1],fits=total,failures,
      failure_rate=p,failure_lower=max(0,center-half),failure_upper=min(1,center+half),
      median_ms=1000*median(group$seconds),p90_ms=1000*unname(quantile(group$seconds,.9)))
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result[order(result$scenario,result$nu,result$method),]
}

# Paired contamination experiment: every optimizer and fixed nu sees the same sample.
# A common t5 baseline is contaminated at nested positions by N(15,1) observations.
# Both data-based and displaced starts are tested; no information SEs are used here.
run_contamination_study <- function(repetitions=100,n=500,nu_values=c(2,5,30),
                                    fractions=c(0,.01,.05),seed=47000,tol=1e-7,maxit=200,
                                    cache_file=NULL,source_dir="project_A",force=FALSE) {
  if (repetitions<2 || n<10 || any(fractions<0 | fractions>=1)) stop("Invalid experiment settings")
  settings <- list(repetitions=repetitions,n=n,nu_values=nu_values,fractions=fractions,
                   seed=seed,tol=tol,maxit=maxit,baseline_nu=5,outlier_mean=15,outlier_sd=1)
  source_files <- sort(list.files(source_dir,pattern="\\.R$",full.names=TRUE))
  signature <- list(settings=settings,source_md5=unname(tools::md5sum(source_files)),R=R.version.string)
  if (!is.null(cache_file) && file.exists(cache_file) && !force) {
    saved <- readRDS(cache_file)
    if (identical(saved$signature,signature)) return(saved)
  }
  methods <- list(EM=fit_em,BFGS=fit_bfgs,Newton=fit_newton)
  invisible(fit_all(simulate_t_data(n,0,1,5,seed=seed),5,tol=tol))
  rows <- list()
  for (replication in seq_len(repetitions)) {
    dataset_seed <- seed+replication
    baseline <- simulate_t_data(n,0,1,5,seed=dataset_seed)
    positions <- sample.int(n)
    outliers <- rnorm(n,mean=15,sd=1)
    for (fraction in fractions) {
      x <- baseline; replacements <- floor(n*fraction)
      if (replacements>0) x[positions[seq_len(replacements)]] <- outliers[seq_len(replacements)]
      initial_scale <- IQR(x)/2
      starts <- list("Data-based"=c(median(x),2*log(initial_scale)),
                     "Displaced"=c(median(x)+5*initial_scale,2*log(initial_scale/2)))
      for (nu in nu_values) for (scenario in names(starts)) {
        order <- ((seq_along(methods)+replication-2) %% length(methods))+1
        for (method in names(methods)[order]) {
          tic <- proc.time()[3]
          fit <- tryCatch(methods[[method]](x,nu,start=starts[[scenario]],tol=tol,maxit=maxit),
                          error=function(e) e)
          seconds <- unname(proc.time()[3]-tic)
          if (inherits(fit,"error")) {
            row <- data.frame(mu=NA_real_,v=NA_real_,loglik=NA_real_,converged=FALSE,
                              score=NA_real_,reason=conditionMessage(fit))
          } else {
            score <- score_error(x,c(fit$mu,log(fit$v)),nu)
            certified <- isTRUE(fit$converged) && is.finite(score) && score<=tol
            row <- data.frame(mu=fit$mu,v=fit$v,loglik=fit$loglik,converged=certified,score=score,
                              reason=if(certified) "Certified" else "Stopping rule / score not certified")
          }
          rows[[length(rows)+1]] <- cbind(data.frame(replication,dataset_seed,fraction,replacements,
            nu,scenario,method,seconds,start_mu=starts[[scenario]][1],start_logv=starts[[scenario]][2]),row)
        }
      }
    }
  }
  result <- list(signature=signature,settings=settings,fits=do.call(rbind,rows),
                 completed_at=format(Sys.time(),tz="Europe/Copenhagen",usetz=TRUE),session=sessionInfo())
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file),recursive=TRUE,showWarnings=FALSE)
    saveRDS(result,cache_file)
  }
  result
}

# Summarize within each nu/contamination/start/method cell. Failures remain in timings.
contamination_summary <- function(study) {
  groups <- split(study$fits,interaction(study$fits[c("nu","fraction","scenario","method")],drop=TRUE))
  rows <- lapply(groups,function(group) {
    total <- nrow(group); failures <- sum(!group$converged); p <- failures/total; z <- qnorm(.975)
    center <- (p+z^2/(2*total))/(1+z^2/total)
    half <- z*sqrt(p*(1-p)/total+z^2/(4*total^2))/(1+z^2/total)
    data.frame(nu=group$nu[1],percent=100*group$fraction[1],scenario=group$scenario[1],method=group$method[1],
      fits=total,failures,failure_rate=p,failure_lower=max(0,center-half),failure_upper=min(1,center+half),
      median_ms=1000*median(group$seconds),p90_ms=1000*unname(quantile(group$seconds,.9)))
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result[order(result$scenario,result$nu,result$percent,result$method),]
}

# Agreement is meaningful only within one dataset, fitted nu and contamination level.
# Certified fits from both starts supply a best-attained (not proven global) reference.
contamination_agreement <- function(study) {
  groups <- split(study$fits,interaction(study$fits[c("replication","nu","fraction")],drop=TRUE))
  rows <- lapply(groups,function(group) {
    certified <- group[group$converged,]
    if (nrow(certified)==0) return(data.frame(nu=group$nu[1],percent=100*group$fraction[1],
      replication=group$replication[1],certified=0,mu_spread=NA_real_,sigma_spread=NA_real_,
      loglik_spread=NA_real_,uncertified_gap=NA_real_,uncertified_score=NA_real_))
    reference <- max(certified$loglik)
    failed <- group[!group$converged & is.finite(group$loglik),]
    data.frame(nu=group$nu[1],percent=100*group$fraction[1],replication=group$replication[1],
      certified=nrow(certified),mu_spread=diff(range(certified$mu)),
      sigma_spread=diff(range(sqrt(certified$v))),loglik_spread=diff(range(certified$loglik)),
      uncertified_gap=if(nrow(failed)>0) max(reference-failed$loglik) else NA_real_,
      uncertified_score=if(nrow(failed)>0) max(failed$score) else NA_real_)
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result[order(result$nu,result$percent,result$replication),]
}
