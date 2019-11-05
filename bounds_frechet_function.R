library(gtools)   # for rdirichlet()
library(MASS)     # for mvrnorm()
library(Matrix)
library(lpSolve)  # for sensitivity analysis
library(reshape2)
library(plyr)

`%.%` <- paste0

#############################
## docs for bounds.wrapper ##
#############################

#' Calculate bounds on average choice-specific treatment effects (ACTE) in a
#' patient preference trial with design described in Knox, Yamamoto, Baum, and
#' Berinsky (n.d.), "Design, Identification, and Sensitivity Analysis for
#' Patient Preference Trials."
#'
#' @param data Data frame with the following columns:, 'S' (stated prefs), 'C'
#'     (actual choices), 'A' \itemize{ \item{\code{D}}{ Treatment arm
#'     assignment, with values 0=free choice or 1=forced choice }
#'     \item{\code{S}}{ Stated preference, with values 1:J } \item{\code{C}}{
#'     Actual choice, with values in 1:J if D=0, or NA if D=1 } \item{\code{A}}{
#'     Assigned treatment, with values in 1:J if D=1, or A=C if D=0 }
#'     \item{\code{Y}}{ Continuous numeric outcome } }
#' @param effects.mat Logical matrix of dimensions JxJ, representing the
#'     treatment effects to compute. Rows represent control conditions and
#'     columns represent treatment conditions. For example, to compute the
#'     conditional average treatment effect Y(3) - Y(2)---i.e., changing from
#'     A=2 to A=3---set \code{effects.mat[2, 3]} to \code{TRUE}. To calculate
#'     treatment effects for all choose(J, 2) possible manipulations, set equal
#'     to \code{matrix(TRUE, J, J)}.
#' @param posterior Logical. Simulate from posterior?
#' @param n_dir_MC Integer. Number of Dirichlet draws for posterior of
#'     stated/actual choices. Together, \code{n_dir_MC * n_norm_MC} determine
#'     the number of Monte Carlo simulations from the posterior.
#' @param n_norm_MC Integer. Number of multivariate normal draws (per Dirichlet
#'     draw) for posterior of bounds, conditional on stated/actual
#'     choices. Together, \code{n_dir_MC * n_norm_MC} determine the number of
#'     Monte Carlo simulations from the posterior.
#' @param quantiles Numeric vector with values in [0, 1]. Quantiles of posterior to report.
#' @param sensitivity Logical. Compute sensitivity analysis?
#' @param rhos Numeric vector. Values of sensitivity parameter to compute.



######################################################
## helper functions for variance of integrated ecdf ##
######################################################

integrate_ecdf = function(F_hat, from, to, knots=NULL){

    if (from==to)
        return(0)

    if (from > to)
        stop("'from' must be <= 'to'")

    if (is.null(NULL)){
        if ("stepfun" %in% class(F_hat)){
            knots_trim = knots(F_hat)
        } else if ("knots" %in% names(attributes(F_hat))){
            knots_trim = attr(F_hat,"knots")
        }
    } else {
        knots_trim = knots
    }

    ## this assumes that 'from' and 'to' are in 'knots'
    knots_trim = knots_trim[knots_trim <= to]
    knots_trim = knots_trim[knots_trim > from]
    knots_trim = c(from, knots_trim)

    if (length(knots_trim)==1){
        return(0)
    }

    dx = diff(knots_trim)
    if ("stepfun" %in% class(F_hat)){
        h_x = F_hat(knots_trim[-length(knots_trim)])
    } else if ("knots" %in% names(attributes(F_hat))){
        h_x = sapply(knots_trim[-length(knots_trim)],F_hat)
    }

    out = sum(h_x * dx)
    if (to > max(knots_trim))
        out = out + (to - max(knots_trim)) * F_hat(max(knots_trim))

    return(out)
}
## checking integrate_ecdf():
## n = 100
## x = runif(n)
## mean(x) == diff(range(x)) - integrate_ecdf(ecdf(x), min(x), max(x)) + min(x)
## should be nearly identical


## This is the inner integral in the calculation
##   of the variance of the bounds
## The integral of the empirical survival
##   function is simply the total area
##   of a series of rectangles with width
##   dx[i] and height h_x[i]
integrate_inner = function(F_hat, from, to, knots=NULL){

    if (from==to)
        return(0)

    if (from > to)
        stop("'from' must be <= 'to'")

    if (is.null(NULL)){
        if ("stepfun" %in% class(F_hat)){
            knots_trim = knots(F_hat)
        } else if ("knots" %in% names(attributes(F_hat))){
            knots_trim = attr(F_hat,"knots")
        }
    } else {
        knots_trim = knots
    }

    knots_trim = knots_trim[knots_trim <= to]
    knots_trim = knots_trim[knots_trim > from]
    knots_trim = c(from, knots_trim)

    dx = diff(knots_trim)
    h_x = 1 - F_hat(knots_trim[-length(knots_trim)])

    return(sum(h_x * dx))
}

## This is the outer integral in the calculation
##   of the variance of the bounds
## Let h_y be the value of the outer integrand,
##   F_hat(y) * integrate_esurv(F_hat,y,Q_b) ## (empirical survival function)
## h_y varies linearly within subregions
##   lying between knots, i.e., between observed
##   values of the outcome variable.
## Then dy is the size of the subregion
## The outer integral is then equivalent
##   to the total area of a sequence of
##   trapezoids with widths dy[i] and
##   linearly varying heights
## h_y_bar is the average value of h_y over
##   a given subregion
## The outer integral is then equivalent to
##   the total area of a sequence of rectangles
##   with widths dy[i] and heights h_y[i]
integrate_outer = function(F_hat, from, to, knots=NULL){
    if (is.null(NULL)){
        if ("stepfun" %in% class(F_hat)){
            knots_trim_outer = knots(F_hat)
        } else if ("knots" %in% names(attributes(F_hat))){
            knots_trim_outer = attr(F_hat,"knots")
        }
    } else {
        knots_trim_outer = knots
    }

    if (from==to)
        return(0)

    if (from > to)
        stop("'from' must be <= 'to'")

    ## integrate covariance of ECDFs with intersecting integral domains
    knots_trim_outer = knots_trim_outer[knots_trim_outer <= to]
    knots_trim_outer = knots_trim_outer[knots_trim_outer > from]
    knots_trim_outer = c(from, knots_trim_outer)

    if (length(knots_trim_outer)==1){
        return(0)
    }

    h_y_bar = sapply(1:(length(knots_trim_outer)-1),function(i){
        F_hat(knots_trim_outer[i]) *
            (integrate_inner(F_hat, knots_trim_outer[i], to, knots_trim_outer[i:length(knots_trim_outer)]) +
             integrate_inner(F_hat, knots_trim_outer[i+1], to, knots_trim_outer[(i+1):length(knots_trim_outer)])) / 2
    })
    dy = diff(knots_trim_outer)

    return(sum(h_y_bar * dy))
}
## checking integrate_outer(), should be identical:
## n = 100
## x = runif(n)
## var(x) == 2 * n / (n-1) * integrate_outer(ecdf(x), 0, 1)  #

## ## More checks. Compare
## ##   (1) the variance of the estimated mean calculated
## ##       via the integral of the survival function, and
## ##   (2) the variance of the estimated mean calculated
## ##       via standard methods
## set.seed(02139)
## n_test = 100
## test = rnorm(n_test)
## f_test = ecdf(test)
## integrate_outer(f_test,
##                 from=min(test),
##                 to=max(test)) * 2 / (n_test-1)
## var(test) / n_test
## ## the two are identical to nine digits, after applying
## ##   n/(n-1) correction to the integral method



integrate_ecdf_cov = function(F_hat, from1, to1, from2, to2, knots){
    if (from2 < from1)
        stop("reverse order of integrals")

    breakpoints = sort(unique(c(from1, to1, from2, to2)))
    regions = lapply(1:(length(breakpoints)-1), function(i) breakpoints[i + 0:1])
    n_region = length(regions)

    regions1 = sapply(regions, function(region){
        from1 <= region[1] & to1 >= region[2]
    })
    regions2 = sapply(regions, function(region){
        from2 <= region[1] & to2 >= region[2]
    })

    running_sum = 0
    for (regionA in which(regions1)){
        for (regionB in which(regions2 & 1:n_region >= regionA)){
            if (regionA==regionB){
                running_sum = running_sum + 2 * integrate_outer(F_hat, from=regions[[regionA]][1], to=regions[[regionA]][2], knots)
            } else {
                running_sum = running_sum +
                    integrate_ecdf(F_hat, from=regions[[regionA]][1], to=regions[[regionA]][2], knots) *
                    (diff(regions[[regionB]]) - integrate_ecdf(F_hat, from=regions[[regionB]][1], to=regions[[regionB]][2], knots))
            }
        }
    }

    return(running_sum)
}



###########################################
## calculate shortest posterior interval ##
##   that fully contains 95% mass        ##
###########################################

hpd_bounds <- function(min_MC, max_MC, alpha){
  sol = tryCatch(
  {
    optim(
      unname(quantile(max_MC, alpha + (1-alpha)/2)), # init
      fn = function (cihi){
        x = min_MC
        x[max_MC > cihi] = min(x)
        cilo = quantile(x, 1 - alpha)
        return(cihi - cilo)
      },
      lower = quantile(max_MC, alpha),
      upper = quantile(max_MC, alpha + (1-alpha)/2),
      method = 'L-BFGS-B'
    )
  },
  error = function(e){
    NULL
  })
  if (!is.null(sol)){
    return(c(
      min_cilo = sol$par - sol$value,
      max_cihi = sol$par
    ))
  } else {
    return(c(
      min_cilo = quantile(min_MC, (1-alpha)/2, na.rm=TRUE),
      max_cihi = quantile(max_MC, alpha + (1-alpha)/2, na.rm=TRUE)
    ))
  }
}



#############################################
## bounds calculation, including posterior ##
#############################################

bounds.frechet.calc = function(F_Y_sa, F_Y_scc, F_diff_Y_sa, p_s.c, p_c.s, n_sa_force, n_sc, y.min, y.max, posterior = FALSE, sensitivity = FALSE){

    J = nrow(p_c.s)

    min_sca = max_sca =
        array(NA,rep(J,3),dimnames=list(S="S." %.% 1:J,C="C." %.% 1:J,A="A." %.% 1:J))

    ## if (sensitivity){
        naive_sa = matrix(NA, J, J, dimnames=list(C="C." %.% 1:J,A="A." %.% 1:J))
    ## }

    if (posterior){
        min_sca.var = max_sca.var = minmax_sca.cov =
            array(NA,rep(J,3),dimnames=list(S="S." %.% 1:J,C="C." %.% 1:J,A="A." %.% 1:J))
        if (sensitivity){
            naive_sa.var = matrix(NA, J, J, dimnames=list(C="C." %.% 1:J,A="A." %.% 1:J))
            naivemin_sca.cov = naivemax_sca.cov =
                array(0,rep(J,3),dimnames=list(S="S." %.% 1:J,C="C." %.% 1:J,A="A." %.% 1:J))
        }
    }

    for (s in 1:J){
        for (c in 1:J){
            for (a in (1:J)[-c]){

                b = p_c.s[s,c] / (1 - p_c.s[s,a])

                ## ignore warnings about "no non-missing to min" (returns Inf);
                ##   this is handled in integrate_ecdf
                Q_b.ind = suppressWarnings(min(which(attr(F_diff_Y_sa[[s,a]],"values") > b)))
                Q_b = attr(F_diff_Y_sa[[s,a]],"knots")[Q_b.ind]

                Q_1mb.ind = suppressWarnings(min(which(attr(F_diff_Y_sa[[s,a]],"values") > 1-b)))
                Q_1mb = attr(F_diff_Y_sa[[s,a]],"knots")[Q_1mb.ind]

                ## When Pr(C=a|S=s) is extremely high, sampling variability can lead to
                ## F_diff_Y_sa going outside [0,1]. If this happens, the approach we use here
                ## to calculate the lower bound no longer guarantees it will be less than the upper bound
                lb = Q_b - integrate_ecdf(F_diff_Y_sa[[s,a]],y.min,Q_b) / b
                ub = Q_1mb + (y.max - Q_1mb) / b - integrate_ecdf(F_diff_Y_sa[[s,a]],Q_1mb,y.max) / b
                min_sca[s,c,a] = min(lb, ub)
                max_sca[s,c,a] = max(lb, ub)

            }

            min_sca[s,c,c] = max_sca[s,c,c] = y.max - integrate_ecdf(F_Y_scc[[s,c]], from=y.min, to=y.max)

            if (posterior & !sensitivity){
                min_sca.var[s,c,c] = max_sca.var[s,c,c] =
                    integrate_outer(F_Y_scc[[s,c]], from=y.min, to=y.max) * 2 / (n_sc[s,c] - 1)
                minmax_sca.cov[s,c,c] = min_sca.var[s,c,c]
            }

        }

        for (a in 1:J){
            naive_sa[s,a] = y.max - integrate_ecdf(F_Y_sa[[s,a]], y.min, y.max)
        }

        if (posterior){

            ## variance-covariance matrix of all parameters given S_i=s, A_i=a
            vcov_sa = matrix(list(), J, J)

            for (s in 1:J){
                for (a in 1:J){
                    varnames = c(as.character(rbind("min_C." %.% 1:J, "max_C." %.% 1:J)), "naive", "point")
                    ## min/max_C.c are bounds on pi_sca; naive is E[Y(a) | S=s]; point is E[Y(a) | S=s, C=a]

                    ## ignore warnings about "no non-missing to min" (returns Inf);
                    ##   this is handled in integrate_ecdf
                    Q_b = sapply(1:J, function(c){
                        b = p_c.s[s,c] / (1 - p_c.s[s,a])
                        Q_b.ind = suppressWarnings(min(which(attr(F_diff_Y_sa[[s,a]],"values") > b)))
                        return(attr(F_diff_Y_sa[[s,a]],"knots")[Q_b.ind])
                    })
                    Q_1mb = sapply(1:J, function(c){
                        b = p_c.s[s,c] / (1 - p_c.s[s,a])
                        Q_1mb.ind = suppressWarnings(min(which(attr(F_diff_Y_sa[[s,a]],"values") > 1-b)))
                        return(attr(F_diff_Y_sa[[s,a]],"knots")[Q_1mb.ind])
                    })

                    ## varname[i] calculated by integrating F_Y_scc and F_Y_sa
                    ##   between from[i] and to[i]
                    from = c(as.numeric(rbind(y.min, Q_1mb)), y.min, y.min)
                    to = c(as.numeric(rbind(Q_b, y.max)), y.max, y.max)

                    coef_force = c(rep(1/p_c.s[s,1:J], each=2), 1, 0)
                    coef_free = c(rep(p_c.s[s,a]/p_c.s[s,1:J], each=2), 0, 1)

                    drop = grep("C\\." %.% a, varnames)  # no bounds--these are point-identified
                    coef_force[drop] = coef_free[drop] = NA

                    coef_force.pair = coef_force %o% coef_force
                    coef_free.pair = coef_free %o% coef_free

                    cov_force = cov_free = matrix(0, 2*(J+1), 2*(J+1))

                    ijs = which(tril(coef_force.pair) > 0, arr.ind=T)
                    for (ij.ind in 1:nrow(ijs)){
                        ij = ijs[ij.ind,]
                        i = ij[1]; j = ij[2]
                        if (from[i] > from[j] | (from[i] == from[j] & to[i] > to[j])){
                            i = ij[2]; j = ij[1]
                        }

                        cov_force[i,j] = cov_force[j,i] = coef_force.pair[i,j] / (n_sa_force[s,a] - 1) *
                            integrate_ecdf_cov(F_Y_sa[[s,a]], from[i], to[i], from[j], to[j], attr(F_diff_Y_sa[[s,a]],"knots"))
                    }

                    ijs = which(tril(coef_free.pair) > 0, arr.ind=T)
                    for (ij.ind in 1:nrow(ijs)){
                        ij = ijs[ij.ind,]
                        i = ij[1]; j = ij[2]
                        if (from[i] > from[j] | (from[i] == from[j] & to[i] > to[j])){
                            i = ij[2]; j = ij[1]
                        }

                        cov_free[i,j] = cov_free[j,i] = coef_free.pair[i,j] / (n_sc[s,a] - 1) *
                            integrate_ecdf_cov(F_Y_scc[[s,a]], from[i], to[i], from[j], to[j], attr(F_diff_Y_sa[[s,a]],"knots"))
                    }

                    vcov_sa[[s,a]] = cov_force + cov_free
                    vcov_sa[[s,a]][is.na(coef_force.pair)] = NA
                    rownames(vcov_sa[[s,a]]) = colnames(vcov_sa[[s,a]]) = varnames

                }
            }
        }
    }

    out <- list()

    out$max_ca = apply(max_sca * p_s.c %o% rep(1,J),c(2,3),sum)
    out$min_ca = apply(min_sca * p_s.c %o% rep(1,J),c(2,3),sum)
    dimnames(out$min_ca) <- dimnames(out$max_ca) <- list(C = NULL, A = NULL)
    out$strata <- melt(out$min_ca, value.name = 'min')
    out$strata$max <- melt(out$max_ca)$value

    effects.colnames <- c('C', 'treat', 'ctrl', 'min', 'max')
    effects.ind = which(effects.mat, arr.ind=T)
    colnames(effects.ind) = c("ctrl", "treat")
    effects.min = effects.max =
        array(NA, dim=c(J, nrow(effects.ind)), dimnames=list(C=1:J, effect=NULL))
    for (effect in 1:nrow(effects.ind)){
        effects.min[,effect] = out$min_ca[,effects.ind[effect, "treat"]] - out$max_ca[,effects.ind[effect, "ctrl"]]
        effects.max[,effect] = out$max_ca[,effects.ind[effect, "treat"]] - out$min_ca[,effects.ind[effect, "ctrl"]]
    }
    effects = melt(effects.min)
    colnames(effects)[colnames(effects)=="value"] = "min"
    effects$max = melt(effects.max)$value
    effects$ctrl = effects.ind[effects$effect,"ctrl"]
    effects$treat = effects.ind[effects$effect,"treat"]
    effects$effect = NULL
    out$effects <- effects[,effects.colnames]

    out$naive <- list()
    dimnames(naive_sa) <- list(C = NULL, A = NULL)
    out$naive$strata <- melt(naive_sa, value.name = 'naive')
    colnames(out$naive$strata)[match('S', colnames(out$naive$strata))] <- 'C'
    out$naive$effects <- out$effects
    out$naive$effects$naive <-
        naive_sa[cbind(out$effects$C, out$effects$treat)] -
        naive_sa[cbind(out$effects$C, out$effects$ctrl)]
    out$naive$effects <- out$naive$effects[,-match(c('min', 'max'), colnames(out$naive$effects))]

    out$naive_sa = naive_sa
    out$min_sca = min_sca
    out$max_sca = max_sca

    if (posterior){
        out$vcov_sa = vcov_sa
    }

    return(out)
}

sens.calc = function(min_ca,
                     max_ca,
                     naive_sa,
                     p_c.s,
                     p_s.c,
                     effects.mat,
                     rhos,
                     verbose = 0){

    J = nrow(p_c.s)

    ## create containers for solutions
    out = list(min_ca = array(NA, dim=c(J, J, length(rhos)),
                              dimnames=list(C=1:J, A=1:J, rho=NULL)),
               max_ca = array(NA, dim=c(J, J, length(rhos)),
                              dimnames=list(C=1:J, A=1:J, rho=NULL))
               )
    if (any(effects.mat)){
        effects.ind = which(effects.mat, arr.ind = TRUE)
        colnames(effects.ind) = c("ctrl", "treat")
        out[['effects.min']] = out[['effects.max']] =
            array(NA, dim=c(J, length(rhos), nrow(effects.ind)),
                  dimnames=list(C=1:J,
                                rho=NULL,
                                effect=NULL))
    }

    for (rho.ind in length(rhos):1){ # sensitivity

        rho <- rhos[rho.ind]
        ## compute strata sensitivity results
        out$min_ca[,,rho.ind] <- pmax(min_ca, naive_sa - rho)
        out$max_ca[,,rho.ind] <- pmin(max_ca, naive_sa + rho)

        for (c in 1:J){
            for (a in (1:J)){
                for (a2 in which(effects.mat[,a])){
                    ## store acte sensitivity results
                    effect.ind = which(effects.ind[,'ctrl'] == a2 &
                                       effects.ind[,'treat'] == a)
                    out$effects.min[,rho.ind, effect.ind] <-
                        out$min_ca[, a, rho.ind] - out$max_ca[, a2, rho.ind]
                    out$effects.max[,rho.ind, effect.ind] <-
                        out$max_ca[, a, rho.ind] - out$min_ca[, a2, rho.ind]
                }
            }
        }
    }

    strata.colnames <- c('C', 'A', 'rho', 'min', 'max')
    effects.colnames <- c('C', 'treat', 'ctrl', 'rho', 'min', 'max')

    out$strata <- melt(out$min_ca, value.name = 'min')
    out$strata$max <- melt(out$max_ca)$value
    out$strata$rho <- rhos[out$strata$rho]
    out$strata <- out$strata[,strata.colnames]

    out$effects <- melt(out$effects.min, value.name = 'min')
    out$effects$max <- melt(out$effects.max)$value
    out$effects$rho <- rhos[out$effects$rho]
    out$effects$treat = effects.ind[out$effects$effect,"treat"]
    out$effects$ctrl = effects.ind[out$effects$effect,"ctrl"]
    out$effects <- out$effects[,effects.colnames]

    ## drop values of rho that are estimated to be inconsistent with data
    out$strata$min[out$strata$min > out$strata$max] <- NA
    out$strata$max[is.na(out$strata$min)] <- NA
    rho.min <- daply(na.omit(out$strata), c('C', 'A'), function(x){
        min(x$rho)
    })

    drop <-
        out$effects$rho < rho.min[cbind(out$effects$C, out$effects$treat)] |
        out$effects$rho < rho.min[cbind(out$effects$C, out$effects$ctrl)]
    out$effects$min[drop] <- NA
    out$effects$max[drop] <- NA

    return(out)

}



#############################
## wrapper bounds function ##
#############################

bounds.frechet = function(data,
                          effects.mat = NULL,
                          posterior = FALSE,
                          n_dir_MC = 100,
                          n_norm_MC = 10,
                          quantiles = c(.025,.975),
                          sensitivity = FALSE,
                          rhos = Inf,
                          hpd = FALSE,
                          alpha = .95,
                          verbose = 0  ## currently has no effect
                          ) {

    if (!Inf %in% rhos){
        rhos <- c(rhos, Inf)
    }
    if (length(rhos) == 1){
        if (sensitivity){
            warning('no finite values for rhos supplied; setting sensitivity = FALSE')
            sensitivity <- FALSE
        }
    } else {
        if (!sensitivity){
            warning('sensitivity set to FALSE; ignoring rhos')
            rhos <- Inf
        }
    }

    y.min = min(data$Y)
    y.max = max(data$Y)

    ## sum(1{S=s})
    n_s = table(S=data$S)

    ## sum(1{S=s,C=c})
    n_sc = table(S=data$S,C=data$C)

    ## sum(1{S=s,A=a,D=1})
    n_sa_force = table(S=data$S[data$D==1],A=data$A[data$D==1])

    ## Pr(S=s)
    p_s = table(data$S) / nrow(data)

    ## Pr(C=c|S=s)
    p_c.s = t(sapply(1:J,function(s){
        out = table(data[data$S==s & data$D==0,"C"]) / length(data[data$S==s & data$D==0,"C"])
        out = out[as.character(1:J)]
        out[is.na(out)] = 0
        return(out)
    }))
    names(dimnames(p_c.s)) = c("S","C")

    ## Pr(S=s,C=c)
    p_sc_free = n_sc / sum(n_sc)
    p_sc_pool = p_s %o% rep(1,J) * p_c.s # pool stated pref information from forced condition
    dimnames(p_sc_pool) = dimnames(p_sc_free)

    ## Pr(S=s|C=c)
    p_s.c = apply(p_sc_pool,2,function(x){
        out = x / sum(x)
    })
    dimnames(p_s.c) = dimnames(p_sc_pool)

    ## observed outcomes

    ## pi_scc[s,c] == E[Y(c) | S=s,C=c,A=c]
    ## y_sa[s,a] == E[Y(a) | S=s,A=a]

    ## from forced condition: E[Y|S=s,A=a,D=1]
    y_sa = t(sapply(1:J,function(s){
        out = sapply(1:J,function(a){
            mean(data$Y[which(data$S==s & data$A==a & data$D==1)],na.rm=T)
        })
        return(out)
    }))
    dimnames(y_sa) = list(S=NULL,A=NULL)

    ## from free condition: E[Y|S=s,C=c,A=c,D=0]
    pi_scc = t(sapply(1:J,function(s){
        out.s = sapply(1:J,function(c){
            out.c = mean(data$Y[which(data$S==s & data$C==c & data$D==0)],na.rm=T)
            out.c[is.na(out.c)] = 0
            return(out.c)
        })
        return(out.s)
    }))
    dimnames(pi_scc) = list(S=NULL,A=NULL)

    ## ## guide to functions:
    ##
    ## F_Y_sa[[s,a]](y)                  Pr(Y(a)<=y | S=s)
    ## F_Y_scc[[s,c]](y)                 Pr(Y(c)<=y | S=s,C=c)
    ## F_diff_Y_sa[[s,c]](y)             Pr(Y(a)<=y | S=s,C!=a)
    ##
    ## Up until this point, we've been dealing with strata-specific parameters:
    ##   for example, size of the subgroup with S=s and C=c
    ##   or mean potential outcome under treatment 'a' among those with S=s.
    ## These are stored in matrix/array format and accessed as
    ##   p_sc_pool[s,c] and y_sa[s,a], respectively.
    ##
    ## Below we require strata-specific functions,
    ##   e.g., conditional ECDFs among various subgroups.
    ## The number of parameters is unknown and may differ between strata.
    ##
    ## To preserve the matrix-like subscript notation,
    ##   we're storing these functions in list matrices.
    ## For example, the observed ECDF 'Pr(Y(a)<=y | S=s)' is stored in 'F_Y_sa';
    ##   it is accessed by F_Y_sa[[s,a]], and its value at 'y' is 'F_Y_sa[[s,a]](y)'
    ##
    ## Later, we will need to integrate these functions.
    ## Some of them (i.e., the ECDFs) are of class 'stepfun', while
    ##   others operate on these functions (e.g. take differences).
    ## Class 'stepfun' has accessor function 'knots()' which makes integration easy,
    ##   but `+.stepfun` and `*.stepfun` are unfortunately not written.
    ## To preserve similar functionality without hand-rolling these methods,
    ##   whenever a child function is generated, we extract this information
    ##   and save it in a `knots` attribute of the child function.

    ## Pr(Y(a)<=y | S=s)
    F_Y_sa = t(sapply(1:J,function(s){
        sapply(1:J,function(a){
            ecdf(data$Y[which(data$S==s & data$A==a & data$D==1)])
        })
    }))
    dimnames(F_Y_sa) = list(S=NULL,A=NULL)

    ## Pr(Y(c)<=y | S=s,C=c)
    F_Y_scc = t(sapply(1:J,function(s){
        sapply(1:J,function(c){
            ecdf(data$Y[which(data$S==s & data$C==c & data$D==0)])
        })
    }))
    dimnames(F_Y_scc) = list(S=NULL,C=NULL)

    ## Pr(Y(a)<=y | S=s, C!=c)
    F_diff_Y_sa =  t(sapply(1:J,function(s){
        sapply(1:J,function(a){
            ## evaluate 's' and 'a' in local context (single iteration of sapply)
            ## Otherwise, they are evaluated in the parent environment, of 'sapply' itself
            ##   meaning that s=J and a=J (ending values in the parent)
            env_sa = new.env()
            assign("s",s,envir=env_sa) # save current value of 's' into env_sa, where it will remain unchanged by future iterations
            assign("a",a,envir=env_sa) # save current value of 'a' into env_sa, where it will remain unchanged by future iterations
            f_sa = function(y){
                out = (F_Y_sa[[s,a]](y) - F_Y_scc[[s,a]](y) * p_c.s[s,a]) / (1 - p_c.s[s,a])
                return(out)
            }
            environment(f_sa) = env_sa
            attr(f_sa,"knots") = sort(unique(c(y.min,y.max,knots(F_Y_scc[[s,a]]),knots(F_Y_sa[[s,a]]))))
            attr(f_sa,"values") = sapply(attr(f_sa,"knots"),f_sa)
            return(f_sa)
        })
    }))
    dimnames(F_diff_Y_sa) = list(S=NULL,A=NULL)

    ## this assumes that the estimated F_diff_Y_sa is monotonic, which is
    ## necessary to estimate the variance of the bounds. Asymptotically this is
    ## true but for finite samples may not hold. In these cases, the inverse
    ## function Q(b) is taken to be the first value at which F_diff_Y_sa(Q(b)) >= b
    ##
    ## then:
    ##   - the integral of G_Y_sca.lo from y.min to y.max is simplified to an
    ##     integral of F_diff_Y_sa from y.min to Q_b
    ##   - the integral of G_Y_sca.hi from y.min to y.max is simplified to an
    ##     integral of F_diff_Y_sa from Q_1mb to y.max
    ##
    ## implicitly, this assumes that
    ##   - G_Y_sca.lo = 0 for all y > Q_b
    ##     (guaranteed, since Q breaks ties by choosing lowest)
    ##   - G_Y_sca.hi = 1 for all y > Q_1mb
    ##     (not guaranteed in finite samples, since G_Y_sca.hi may not be nondecreasing)

    ## point estimates
    bounds = bounds.frechet.calc(F_Y_sa,
                                 F_Y_scc,
                                 F_diff_Y_sa,
                                 p_s.c,
                                 p_c.s,
                                 n_sa_force,
                                 n_sc,
                                 y.min,
                                 y.max,
                                 posterior = posterior,
                                 sensitivity)

    if (sensitivity){
        sens = sens.calc(bounds$min_ca,
                         bounds$max_ca,
                         bounds$naive_sa,
                         p_c.s,
                         p_s.c,
                         effects.mat,
                         rhos)

        ## if estimated bounds do not exist for some values of rho,
        ##   drop posterior draws for those rhos too
        drop.strata <- which(is.na(as.matrix(sens$strata)), arr.ind = TRUE)
        drop.effects <- which(is.na(as.matrix(sens$effects)), arr.ind = TRUE)

    }

    ## posterior
    if (posterior){

        p_c.s_MC = p_s.c_MC = array(NA, c(n_dir_MC*n_norm_MC,rep(J,2)), dimnames=list(draw=NULL,S="S." %.% 1:J,C="C." %.% 1:J))
        min_sca_MC = max_sca_MC =
            array(NA, c(n_dir_MC*n_norm_MC,rep(J,3)), dimnames=list(draw=NULL,S="S." %.% 1:J,C="C." %.% 1:J,A="A." %.% 1:J))
        min_ca_MC = max_ca_MC =
            array(NA, c(n_dir_MC*n_norm_MC,rep(J,2)), dimnames=list(draw=NULL,C="C." %.% 1:J,A="A." %.% 1:J))
        naive_sa_MC =
            array(NA, c(n_dir_MC*n_norm_MC,rep(J,2)), dimnames=list(draw=NULL,S="S." %.% 1:J,A="A." %.% 1:J))

        if (sensitivity){
            sens_strata_MC =
                array(NA,
                      c(n_dir_MC*n_norm_MC, dim(sens$strata)),
                      dimnames=list(draw=NULL,
                                    sens=NULL,
                                    variable=colnames(sens$strata))
                      )
            sens_effects_MC =
                array(NA,
                      c(n_dir_MC*n_norm_MC, dim(sens$effects)),
                      dimnames=list(draw=NULL,
                                    sens=NULL,
                                    variable=colnames(sens$effects))
                      )
        }

        prev_message_nchar = 0
        for (i in 1:n_dir_MC){

            message_to_display = "dirichlet draw: " %.% i %.% " of " %.% n_dir_MC %.% " (" %.% n_norm_MC %.% " normal draws each)"
            cat(rep("\b", prev_message_nchar), message_to_display, sep="")
            prev_message_nchar = nchar(message_to_display)

            ## draw from posterior of Pr(S=s)
            p_s_star = as.numeric(rdirichlet(1,table(data$S)))

            ## draw from posterior of Pr(C=c|S=s)
            p_c.s_star = t(sapply(1:J,function(s){
                rdirichlet(1,table(data$C[data$S==s]))
            }))
            p_c.s_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),,] = rep(1, n_norm_MC) %o% p_c.s_star

            ## converting to Pr(S=s,C=c)
            p_sc_star = p_s_star %o% rep(1,J) * p_c.s_star # pool stated pref information from forced condition
            dimnames(p_sc_star) = dimnames(p_sc_pool)

            ## Pr(S=s|C=c)
            p_s.c_star = apply(p_sc_star,2,function(x){
                out = x / sum(x)
            })
            p_s.c_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),,] = rep(1, n_norm_MC) %o% p_s.c_star

            bounds_MC = bounds.frechet.calc(F_Y_sa,
                                            F_Y_scc,
                                            F_diff_Y_sa,
                                            p_s.c_star,
                                            p_c.s_star,
                                            n_sa_force,
                                            n_sc,
                                            y.min,
                                            y.max,
                                            posterior = TRUE,
                                            sensitivity)

            for (s in 1:J){
                for (a in 1:J){

                    means = c(
                        as.numeric(rbind(
                            bounds_MC$min_sca[s,,a],
                            bounds_MC$max_sca[s,,a]
                        )),
                        bounds_MC$naive_sa[s,a],
                        bounds_MC$max_sca[s,a,a]
                    )
                    vcov = bounds_MC$vcov_sa[[s,a]]
                    vcov[is.na(vcov)] = 0

                    bounds_sca_MC = mvrnorm(n = n_norm_MC,
                                            mu = means,
                                            Sigma = vcov,
                                            tol = 1)

                    min_sca_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),s,-a,a] = bounds_sca_MC[,(1:J)[-a]*2-1]
                    max_sca_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),s,-a,a] = bounds_sca_MC[,(1:J)[-a]*2]
                    min_sca_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),s,a,a] =
                        max_sca_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),s,a,a] =
                        bounds_sca_MC[,2*J+2]

                    naive_sa_MC[((i-1)*n_norm_MC+1):(i*n_norm_MC),s,a] = bounds_sca_MC[,2*J+1]


                }
            }
        }
        cat("\n")

        ## confidence intervals for bounds on strata means

        max_ca_MC = apply(max_sca_MC * p_s.c_MC %o% rep(1,J),c(1,3,4),sum)
        min_ca_MC = apply(min_sca_MC * p_s.c_MC %o% rep(1,J),c(1,3,4),sum)

        strata.min.quantiles = sapply(quantiles, function(q){
            apply(min_ca_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(strata.min.quantiles) <- list(
            C = NULL,
            A = NULL,
            quantile = sprintf("quantile_%0.3f",quantiles)
        )

        strata.max.quantiles = sapply(quantiles, function(q){
            apply(max_ca_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(strata.max.quantiles) <- list(
            C = NULL,
            A = NULL,
            quantile = sprintf("quantile_%0.3f",quantiles)
        )

        strata.quantiles = sapply(quantiles, function(q){
            out = melt(strata.min.quantiles[,,sprintf("quantile_%0.3f", q),drop=FALSE],
                       value.name = 'min')
            out$max = melt(strata.max.quantiles[,,sprintf("quantile_%0.3f", q)])$value
            out$quantile = NULL
            return(as.matrix(out))
        }, simplify='array')
      dimnames(strata.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

      if (hpd){
        strata.ci <- ddply(
          bounds$strata,
          c('A', 'C'),
          function(comb){
            hpd_bounds(min_ca_MC[,comb$C,comb$A],
                       max_ca_MC[,comb$C,comb$A],
                       alpha
                       )
          })
        strata.ci[c('C', 'A', 'min_cilo', 'max_cihi')]
      }

        ## confidence intervals for bounds on treatment effects

        effects.ind = which(effects.mat, arr.ind=T)
        colnames(effects.ind) = c("ctrl", "treat")
        effects.min_MC = effects.max_MC =
            array(NA, dim=c(n_norm_MC*n_dir_MC, J, nrow(effects.ind)),
                  dimnames=list(draw=NULL,
                                C=1:J,
                                effect=NULL))
        for (effect in 1:nrow(effects.ind)){
            effects.min_MC[,,effect] = min_ca_MC[,,effects.ind[effect, "treat"]] - max_ca_MC[,,effects.ind[effect, "ctrl"]]
            effects.max_MC[,,effect] = max_ca_MC[,,effects.ind[effect, "treat"]] - min_ca_MC[,,effects.ind[effect, "ctrl"]]
        }

        effects.min.quantiles = sapply(quantiles, function(q){
            apply(effects.min_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(effects.min.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)
        names(dimnames(effects.min.quantiles))[3] = 'quantile'

        effects.max.quantiles = sapply(quantiles, function(q){
            apply(effects.max_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(effects.max.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)
        names(dimnames(effects.max.quantiles))[3] = 'quantile'

        effects.colnames <- c('C', 'treat', 'ctrl', 'min', 'max')
        effects.quantiles = sapply(quantiles, function(q){
            ## drop=FALSE to preserve the 'effect' index even when only
            ##   one effect (treatment-control pair) is specified
            ## effect index is used to remember treat/ctrl values then dropped
            out = melt(effects.min.quantiles[,,sprintf("quantile_%0.3f", q),drop=FALSE])
            colnames(out)[colnames(out)=="value"] = "min"
            out$max = melt(effects.max.quantiles[,,sprintf("quantile_%0.3f", q)])$value
            out$ctrl = effects.ind[out$effect,"ctrl"]
            out$treat = effects.ind[out$effect,"treat"]
            out$effect = NULL
            out$quantile = NULL
            return(as.matrix(out))
        }, simplify='array')
        dimnames(effects.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)
      effects.quantiles <- effects.quantiles[,effects.colnames,]

      if (hpd){
        effects.ci <- ddply(
          expand.grid(C = 1:J, effect = 1:nrow(effects.ind)),
          c('effect', 'C'),
          function(comb){
            hpd_bounds(effects.min_MC[,comb$C,comb$effect],
                       effects.max_MC[,comb$C,comb$effect],
                       alpha
                       )
          })
        effects.ci$treat <- effects.ind[effects.ci$effect, 'treat']
        effects.ci$ctrl <- effects.ind[effects.ci$effect, 'ctrl']
        effects.ci$effect <- NULL
        effects.ci <- effects.ci[,c('C', 'treat', 'ctrl', 'min_cilo', 'max_cihi')]
      }

        ## confidence intervals for naive estimates of strata means

        strata.naive.quantiles = sapply(quantiles, function(q){
            apply(naive_sa_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(strata.naive.quantiles) <- list(
            C = NULL,
            A = NULL,
            quantile = sprintf("quantile_%0.3f",quantiles)
        )

        strata.naive.quantiles = sapply(quantiles, function(q){
            out = melt(strata.naive.quantiles[,,sprintf("quantile_%0.3f", q),drop=FALSE],
                       value.name = 'naive')
            out$quantile = NULL
            return(as.matrix(out))
        }, simplify='array')
      dimnames(strata.naive.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

      if (hpd){
        strata.naive.ci <- ddply(
          bounds$naive$strata,
          c('A', 'C'),
          function(comb){
            hpd_bounds(naive_sa_MC[,comb$C,comb$A],
                       naive_sa_MC[,comb$C,comb$A],
                       alpha
                       )
          })
        strata.naive.ci <- strata.naive.ci[,c('C', 'A', 'min_cilo', 'max_cihi')]
      }

        ## confidence intervals for naive estimates of treatment effects

        effects.naive_MC =
            array(NA, dim=c(n_norm_MC*n_dir_MC, J, nrow(effects.ind)),
                  dimnames=list(draw=NULL,
                                C=1:J,
                                effect=NULL))
        for (effect in 1:nrow(effects.ind)){
            effects.naive_MC[,,effect] = naive_sa_MC[,,effects.ind[effect, "treat"]] - naive_sa_MC[,,effects.ind[effect, "ctrl"]]
        }

        effects.naive.quantiles = sapply(quantiles, function(q){
            apply(effects.naive_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(effects.naive.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)
        names(dimnames(effects.naive.quantiles))[3] = 'quantile'

        effects.colnames <- c('C', 'treat', 'ctrl', 'naive')
        effects.naive.quantiles = sapply(quantiles, function(q){
            ## drop=FALSE to preserve the 'effect' index even when only
            ##   one effect (treatment-control pair) is specified
            ## effect index is used to remember treat/ctrl values then dropped
            out = melt(effects.naive.quantiles[,,sprintf("quantile_%0.3f", q),drop=FALSE],
                       value.name = 'naive')
            out$ctrl = effects.ind[out$effect,"ctrl"]
            out$treat = effects.ind[out$effect,"treat"]
            out$effect = NULL
            out$quantile = NULL
            return(as.matrix(out))
        }, simplify='array')
        dimnames(effects.naive.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)
      effects.naive.quantiles <- effects.naive.quantiles[,effects.colnames,]

      if (hpd){
        effects.naive.ci <- ddply(
          expand.grid(C = 1:J, effect = 1:nrow(effects.ind)),
          c('effect', 'C'),
          function(comb){
            hpd_bounds(effects.naive_MC[,comb$C,comb$effect],
                       effects.naive_MC[,comb$C,comb$effect],
                       alpha
                       )
          })
        effects.naive.ci$treat <- effects.ind[effects.naive.ci$effect, 'treat']
        effects.naive.ci$ctrl <- effects.ind[effects.naive.ci$effect, 'ctrl']
        effects.naive.ci$effect <- NULL
        effects.naive.ci <- effects.naive.ci[,c('C', 'treat', 'ctrl', 'min_cilo', 'max_cihi')]
      }

        if (sensitivity){

            prev_message_nchar = 0
            sens_ndraws = n_norm_MC*n_dir_MC
            for (draw in 1:sens_ndraws){
                ## progress reporting
                message_to_display = "sensitivity calculations: " %.% draw %.% " of " %.% sens_ndraws
                cat(rep("\b", prev_message_nchar), message_to_display, sep="")
                prev_message_nchar = nchar(message_to_display)
                ind_dir_MC = draw %% n_norm_MC
                if (ind_dir_MC==0){
                    ind_dir_MC = 10
                }

                sens_MC = sens.calc(min_ca_MC[draw,,],
                                    max_ca_MC[draw,,],
                                    naive_sa_MC[draw,,],
                                    p_c.s_MC[draw,,],
                                    p_s.c_MC[draw,,],
                                    effects.mat,
                                    rhos)

                ## if estimated bounds do not exist for some values of rho,
                ##   drop posterior draws for those rhos too
                sens_MC$strata <- as.matrix(sens_MC$strata)
                sens_MC$strata[drop.strata] <- NA
                sens_strata_MC[draw,,] = sens_MC$strata

                sens_MC$effects <- as.matrix(sens_MC$effects)
                sens_MC$effects[drop.effects] <- NA
                sens_effects_MC[draw,,] = sens_MC$effects

            }
            cat("\n")

            sens_strata.quantiles = sapply(quantiles, function(q){
                apply(sens_strata_MC,c(2,3),quantile,q,na.rm=T)
            },simplify="array")
            dimnames(sens_strata.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

          if (hpd){
            sens_strata.ci <- ddply(
              sens$strata,
              c('rho', 'A', 'C'),
              function(comb){
                ind <-
                  sens_strata_MC[1,,'C'] == comb$C &
                  sens_strata_MC[1,,'A'] == comb$A &
                  sens_strata_MC[1,,'rho'] == comb$rho
                min_MC = sens_strata_MC[, ind, 'min']
                max_MC = sens_strata_MC[, ind, 'max']
                drop = is.na(min_MC) | is.na(max_MC)
                if (all(drop)){
                    return(c(min_cilo = NA_real_,
                             max_cihi = NA_real_,
                             failed_lp = mean(drop)
                             ))
                } else {
                  c(hpd_bounds(min_MC[!drop],
                             max_MC[!drop],
                             alpha
                             ),
                    failed_lp = mean(drop)
                    )
                  }
              })
            sens_strata.ci = sens_strata.ci[,c('C', 'A', 'rho', 'min_cilo', 'max_cihi', 'failed_lp')]
          }

          sens_effects.quantiles = sapply(quantiles, function(q){
            apply(sens_effects_MC,c(2,3),quantile,q,na.rm=T)
          },simplify="array")
          dimnames(sens_effects.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

          if (hpd){
            sens_effects.ci <- as.data.frame(t(sapply(
                1:nrow(sens$effects),
                function(i){
                  comb <- sens$effects[i,]
                  ind <-
                    sens_effects_MC[1,,'C'] == comb$C &
                    sens_effects_MC[1,,'treat'] == comb$treat &
                    sens_effects_MC[1,,'ctrl'] == comb$ctrl &
                    sens_effects_MC[1,,'rho'] == comb$rho
                  min_MC = sens_effects_MC[, ind, 'min']
                  max_MC = sens_effects_MC[, ind, 'max']
                  drop = is.na(min_MC) | is.na(max_MC)
                  if (all(drop)){
                    return(cbind(
                      t(as.numeric(comb[c('C', 'treat', 'ctrl', 'rho')])),
                      t(c(min_cilo = NA_real_, max_cihi = NA_real_)),
                      failed_lp = mean(drop)
                    ))
                  } else {
                    cbind(
                      t(as.numeric(comb[c('C', 'treat', 'ctrl', 'rho')])),
                      t(hpd_bounds(min_MC[!drop],
                                   max_MC[!drop],
                                   alpha
                                   )),
                      failed_lp = mean(drop)
                    )
                  }
                }
              )))
            colnames(sens_effects.ci) <- c('C', 'treat', 'ctrl', 'rho', 'min_cilo', 'max_cihi', 'failed_lp')
          }

        }
    }

    rownames(bounds$min_ca) = colnames(bounds$min_ca) = NULL
    rownames(bounds$max_ca) = colnames(bounds$max_ca) = NULL
    names(dimnames(y_sa)) = c("C", "A")

    out = list(strata = bounds$strata,
               effects = bounds$effects,
               naive = bounds$naive,
               n = nrow(data),
               n_sc_free = n_sc,
               n_sa_force = n_sa_force
               )

    if (posterior){
        out$posterior <- list()
        out$posterior$strata.quantiles = strata.quantiles
        out$posterior$effects.quantiles = effects.quantiles
        out$posterior$naive_strata.quantiles = strata.naive.quantiles
        out$posterior$naive_effects.quantiles = effects.naive.quantiles
        if (hpd){
          out$posterior$strata.ci = strata.ci
          out$posterior$effects.ci = effects.ci
          out$posterior$naive_strata.ci = strata.naive.ci
          out$posterior$naive_effects.ci = effects.naive.ci
        }
    }

    if (sensitivity){
        out$sensitivity = sens[c('strata', 'effects')]
        if (posterior){
            out$posterior$sens_strata.quantiles = sens_strata.quantiles
            out$posterior$sens_effects.quantiles = sens_effects.quantiles
            if (hpd){
              out$posterior$sens_strata.ci = sens_strata.ci
              out$posterior$sens_effects.ci = sens_effects.ci
            }
        }
    }

    return(out)
}

