##########################################
# bounds calculation, including variance #
##########################################

bounds.binary.calc = function(p_s, p_c.s, p_sc_pool,
                              y_sa, pi_scc,
                              effects.mat = NULL,
                              sensitivity=FALSE, rhos = Inf,
                              verbose = 0){

    J = length(p_s)
    outcome.levels = 2
    p_c = colSums(p_sc_pool)

    ## indices for observed strata
    sc.comb = as.matrix(expand.grid(c=1:J,s=1:J)[,2:1])
    sa.comb = sc.comb; colnames(sa.comb) = c("s","a")
    ca.comb = sc.comb; colnames(ca.comb) = c("c","a")
    ca.comb = ca.comb[ca.comb[,"c"]!=ca.comb[,"a"],] # C=A not needed for sensitivity
    sca.comb = as.matrix(expand.grid(a=1:J,c=1:J,s=1:J)[,3:1])

    ## indices for principal strata
    rsc.comb = as.matrix(do.call(expand.grid,
                                 c(list(1:J),list(1:J),lapply(1:J,function(j) 0:(outcome.levels-1)))
                                 )[,(J+2):1])
    colnames(rsc.comb) = c("y1","y2","y3","s","c")
    phi.names = "phi_" %.% apply(rsc.comb,1,paste,collapse="")
    phi.num = outcome.levels^J * J^2

    ## containers for solutions
    out = list(min_ca = array(NA, dim=c(J, J, length(rhos)),
                              dimnames=list(C=1:J, A=1:J, rho=NULL)),
               max_ca = array(NA, dim=c(J, J, length(rhos)),
                              dimnames=list(C=1:J, A=1:J, rho=NULL)),
               naive_sa = y_sa
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



    ## define constraints ##

    ## mechanical constraints
    const_sum.mat = rep(1,phi.num)
    const_sum.rhs = 1
    const_sum.dir = "=="

    ## choice probabilities
    const_p_sc.mat = t(sapply(apply(sc.comb,1,paste,collapse=""),function(sc){
        out = rep(0,phi.num)
        out[grep("^phi_..." %.% sc %.% "$",phi.names)] = 1
        return(out)
    }))
    const_p_sc.dir = rep("==",J^2)
    const_p_sc.rhs = p_sc_pool[sc.comb]


    ## observed outcomes from free choice condition
    const_y_free.mat = t(apply(sc.comb,1,function(row){
        s = row["s"]
        c = row["c"]
        regex = c(rep(".",J),s,c)
        regex[c] = 1
        regex = "^phi_" %.% paste(regex,collapse="") %.% "$"
        out = rep(0,phi.num)
        out[grep(regex,phi.names)] = 1
        return(out)
    }))
    const_y_free.dir = rep("==",J^2)
    const_y_free.rhs = pi_scc[sc.comb] * p_sc_pool[sc.comb]


    ## observed outcomes from forced choice condition
    const_y_force.mat = t(apply(sa.comb,1,function(row){
        s = row["s"]
        a = row["a"]
        regex = c(rep(".",J),s,".")
        regex[a] = 1
        regex = "^phi_" %.% paste(regex,collapse="") %.% "$"
        out = rep(0,phi.num)
        out[grep(regex,phi.names)] = 1
        return(out)
    }))
    const_y_force.dir = rep("==",J^2)
    const_y_force.rhs = y_sa[sa.comb] * rep(p_s,each=J)

    ## assumed sensitivity constraints
    if (sensitivity){
        const_sens.mat = t(apply(ca.comb,1,function(row){
            c = row["c"]
            a = row["a"]
            regex = c(rep(".",J),".",c)
            regex[a] = 1
            regex = "^phi_" %.% paste(regex,collapse="") %.% "$"
            out = rep(0,phi.num)
            out[grep(regex,phi.names)] = 1
            return(out)
        }))
        const_sens_min.dir = rep(">=",J^2-J)
        const_sens_max.dir = rep("<=",J^2-J)
        const_sens_base.rhs = y_sa[ca.comb] * p_c[ca.comb[,'c']]
    }


    ## principal strata names
    colnames(const_p_sc.mat) =
        colnames(const_y_free.mat) =
        colnames(const_y_force.mat) = phi.names



    ## define objective function
    for (rho.ind in length(rhos):1){ # sensitivity
        for (obj.c in 1:J){
            for (obj.a in (1:J)){
                for (obj.a2 in c(NA, which(effects.mat[,obj.a]))){

                    rho <- rhos[rho.ind]
                    ## only apply sensitivity constraints for choices/treatments
                    ## of interest
                    sens.ind <- which(ca.comb[,'c'] == obj.c &
                                      ca.comb[,'a'] %in% c(obj.a, obj.a2)
                                      )

                    ## when a2 is NA, bound E[Y(a) | C=c]
                    obj = rep(0,phi.num)
                    regex = c(rep(".",J),".",obj.c)
                    regex[obj.a] = 1
                    regex = "^phi_" %.% paste(regex,collapse="") %.% "$"
                    obj[grep(regex,phi.names)] = 1

                    ## otherwise bound E[Y(a) - Y(a2) | C=c]
                    ## if calculating a treatment effect, subtract
                    ## mean under control
                    if (!is.na(obj.a2)){
                        regex = c(rep(".",J),".",obj.c)
                        regex[obj.a2] = 1
                        regex = "^phi_" %.% paste(regex,collapse="") %.% "$"
                        obj[grep(regex,phi.names)] =
                            obj[grep(regex,phi.names)] - 1
                        if (verbose >= 2){
                            message(sprintf('bounding E[Y(%d)-Y(%d)|C=%d] with rho=%s', obj.a, obj.a2, obj.c, rho))
                        }
                    } else {
                        if (verbose >= 2){
                            message(sprintf('bounding E[Y(%d)|C=%d] with rho=%s', obj.a, obj.c, rho))
                        }
                    }

                    for (direction in c('min', 'max')){
                        ## ## debug
                        ## if (verbose > 5000){
                        ##     cat('  ' %.% direction %.% 'imizing\n')
                        ## }
                        ## ## /debug

                        sol = lp(direction = direction,
                                 objective.in = obj,
                                 const.mat = rbind(const_sum.mat,
                                                   const_p_sc.mat,
                                                   const_y_free.mat,
                                                   const_y_force.mat,
                                                   if (rho < Inf){
                                                       rbind(
                                                           const_sens.mat[sens.ind,],
                                                           const_sens.mat[sens.ind,]
                                                       )
                                                   }
                                                   ),
                                 const.dir = c(const_sum.dir,
                                               const_p_sc.dir,
                                               const_y_free.dir,
                                               const_y_force.dir,
                                               if (rho < Inf){
                                                   c(
                                                       const_sens_min.dir[sens.ind],
                                                       const_sens_max.dir[sens.ind]
                                                   )
                                               }
                                               ),
                                 const.rhs = c(const_sum.rhs,
                                               const_p_sc.rhs,
                                               const_y_free.rhs,
                                               const_y_force.rhs,
                                               if (rho < Inf){
                                                   c(
                                                       const_sens_base.rhs[sens.ind] - rho * p_c[ca.comb[sens.ind, 'c']],
                                                       const_sens_base.rhs[sens.ind] + rho * p_c[ca.comb[sens.ind, 'c']]
                                                   )
                                               }
                                               )
                                 )

                        if (sol$status == 0){
                            if (is.na(obj.a2)){
                                out[[direction %.% '_ca']][obj.c, obj.a, rho.ind] = sol$objval / p_c[obj.c]
                            } else {
                                effect.ind = which(effects.ind[,'ctrl'] == obj.a2 &
                                                   effects.ind[,'treat'] == obj.a)
                                out[['effects.' %.% direction]][obj.c, rho.ind, effect.ind] =
                                    sol$objval / p_c[obj.c]
                            }
                        }
                        if (verbose >= 3){
                            message("  " %.% direction %.% "imization " %.% ifelse(sol$status==0,"succeeded","failed"))
                        }
                    }
                    rm(obj)
                }
            }
        }
    }

    strata.colnames <- c('C', 'A', 'rho', 'min', 'max')
    effects.colnames <- c('C', 'treat', 'ctrl', 'rho', 'min', 'max')

    out$sens <- list()

    out$sens$strata <- melt(out$min_ca, value.name = 'min')
    out$sens$strata$max <- melt(out$max_ca)$value
    out$sens$strata$rho <- rhos[out$sens$strata$rho]
    out$sens$strata <- out$sens$strata[,strata.colnames]

    out$sens$effects <- melt(out$effects.min, value.name = 'min')
    out$sens$effects$max <- melt(out$effects.max)$value
    out$sens$effects$rho <- rhos[out$sens$effects$rho]
    out$sens$effects$treat = effects.ind[out$sens$effects$effect,"treat"]
    out$sens$effects$ctrl = effects.ind[out$sens$effects$effect,"ctrl"]
    out$sens$effects <- out$sens$effects[,effects.colnames]

    out$strata <- out$sens$strata[out$sens$strata$rho == Inf,
                                  -match('rho', strata.colnames)]
    out$effects <- out$sens$effects[out$sens$effects$rho == Inf,
                                    -match('rho', effects.colnames)]
    rownames(out$strata) <- rownames(out$effects) <- NULL

    out$naive <- list()
    out$naive$strata <- melt(y_sa, value.name = 'naive')
    colnames(out$naive$strata)[match('S', colnames(out$naive$strata))] <- 'C'
    out$naive$effects <- out$effects
    out$naive$effects$naive <-
        y_sa[cbind(out$effects$C, out$effects$treat)] -
        y_sa[cbind(out$effects$C, out$effects$ctrl)]
    out$naive$effects <- out$naive$effects[,-match(c('min', 'max'), colnames(out$naive$effects))]


    if (!sensitivity){
        out$sens <- NULL
    } else {
        names(out)[match('sens', names(out))] <- 'sensitivity'
    }

    return(out[c('strata', 'effects', 'naive', 'sensitivity')])

}


#############################
## wrapper bounds function ##
#############################

bounds.binary = function(data,
                         effects.mat = NULL,
                         posterior = FALSE,
                         n_MC = 1000,
                         quantiles = c(.025,.975),
                         sensitivity = FALSE,
                         rhos = Inf,
                         verbose = 0
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

    J = length(na.omit(unique(c(
        data$S, data$C, data$A
    ))))

    if (!all(as.numeric(data$Y) %in% 0:1)){
        stop('outcome levels >3 not implemented')
    }

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
    p_sc_pool = p_s %o% rep(1,J) * p_c.s ## pool stated pref information from forced condition
    dimnames(p_sc_pool) = dimnames(p_sc_free)

    ## Pr(S=s|C=c)
    p_s.c = apply(p_sc_pool,2,function(x){
        out = x / sum(x)
    })
    dimnames(p_s.c) = dimnames(p_sc_pool)



    ## observed outcomes ##
    ## pi_scc[s,c]                     E[Y(c) | S=s,C=c,A=c]
    ## y_sa[s,a]                       E[Y(a) | S=s,A=a]

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

    ## calculate bounds
    bounds = bounds.binary.calc(p_s, p_c.s, p_sc_pool, # preferences
                                y_sa, pi_scc,          # observed outcomes
                                effects.mat,           # effects to calculate
                                sensitivity, rhos,     # sensitivity parameters
                                verbose
                                )

    ## if estimated bounds do not exist for some values of rho,
    ##   drop posterior draws for those rhos too
    drop.strata <- which(is.na(as.matrix(bounds$sensitivity$strata)), arr.ind = TRUE)
    drop.effects <- which(is.na(as.matrix(bounds$sensitivity$effects)), arr.ind = TRUE)

    effects.ind = which(effects.mat, arr.ind = TRUE)
    colnames(effects.ind) = c("ctrl", "treat")

    effects = as.matrix(bounds$effects)

    ## simulate from posterior
    if (posterior){

        strata_MC = array(NA,
                           dim = c(n_MC, dim(bounds$strata)),
                           dimnames = c(draw = list(NULL), dimnames(bounds$strata)))

        effects_MC = array(NA,
                           dim = c(n_MC, dim(bounds$effects)),
                           dimnames = c(draw = list(NULL), dimnames(bounds$effects)))

        naive_strata_MC =
            array(NA, c(n_MC, dim(bounds$naive$strata)),
                  dimnames=c(draw = list(NULL),
                             strata = list(NULL),
                             list(colnames(bounds$naive$strata))
                             )
                  )

        naive_effects_MC =
            array(NA, c(n_MC, dim(bounds$naive$effects)),
                  dimnames=c(draw = list(NULL),
                             effect = list(NULL),
                             list(colnames(bounds$naive$effects))
                             )
                  )

        if (sensitivity){
            sens_strata_MC =
                array(NA, c(n_MC, dim(bounds$sensitivity$strata)),
                      dimnames=c(draw = list(NULL),
                                 strata = list(NULL),
                                 list(colnames(bounds$sensitivity$strata))
                                 )
                      )
            sens_effects_MC =
                array(NA, c(n_MC, dim(bounds$sensitivity$effects)),
                      dimnames=c(draw = list(NULL),
                                 effect = list(NULL),
                                 list(colnames(bounds$sensitivity$effects))
                                 )
                      )
        }

        prev_message_nchar = 0
        for (i in 1:n_MC){

            if (verbose >= 1){
                message_to_display = "posterior draw: " %.% i %.% " of " %.% n_MC
                cat(rep("\b", prev_message_nchar), message_to_display, sep="")
                prev_message_nchar = nchar(message_to_display)
            }

            ## draw from posterior of Pr(S=s)
            p_s_star = as.numeric(rdirichlet(1, n_s))

            ## draw from posterior of Pr(C=c|S=s)
            p_c.s_star = t(sapply(1:J,function(s){
                rdirichlet(1, n_sc[s,])
            }))
            ## p_c.s_MC[i,,] = p_c.s_star

            ## converting to Pr(S=s,C=c)
            p_sc_star = p_s_star %o% rep(1,J) * p_c.s_star # pool stated pref information from forced condition
            dimnames(p_sc_star) = dimnames(p_sc_pool)

            ## Pr(S=s|C=c)
            p_s.c_star = apply(p_sc_star, 2, function(x){
                out = x / sum(x)
            })
            ## p_s.c_MC[i,,] = p_s.c_star

            ## E[Y(a)|S=s]
            y_sa_star = matrix(
                rbeta(J^2,
                      y_sa * n_sa_force,
                      (1 - y_sa) * n_sa_force),
                nrow = J, ncol = J)
            dimnames(y_sa_star) <- dimnames(y_sa)

            ## E[Y(a)|S=s]
            pi_scc_star = matrix(
                rbeta(J^2,
                      pi_scc * n_sc,
                      (1 - pi_scc) * n_sc),
                nrow = J, ncol = J)
            dimnames(pi_scc_star) <- dimnames(pi_scc)

            bounds_star = bounds.binary.calc(p_s_star, p_c.s_star, p_sc_star,
                                             y_sa_star, pi_scc_star,
                                             effects.mat,
                                             sensitivity, rhos,
                                             verbose - 2)

            strata_MC[i,,] <- as.matrix(bounds_star$strata)
            effects_MC[i,,] <- as.matrix(bounds_star$effects)

            if (sensitivity){

                ## if estimated bounds do not exist for some values of rho,
                ##   drop posterior draws for those rhos too
                bounds_star$sensitivity$strata <- as.matrix(bounds_star$sensitivity$strata)
                bounds_star$sensitivity$strata[drop.strata] <- NA
                sens_strata_MC[i,,] = bounds_star$sensitivity$strata

                bounds_star$sensitivity$effects <- as.matrix(bounds_star$sensitivity$effects)
                bounds_star$sensitivity$effects[drop.effects] <- NA
                sens_effects_MC[i,,] = bounds_star$sensitivity$effects

            }

            naive_strata_MC[i,,] = as.matrix(melt(y_sa_star))
            naive_effects = bounds$effects
            naive_effects$naive <-
                y_sa_star[cbind(bounds$effects$C, bounds$effects$treat)] -
                y_sa_star[cbind(bounds$effects$C, bounds$effects$ctrl)]
            naive_effects <- naive_effects[,-match(c('min', 'max'), colnames(naive_effects))]
            naive_effects_MC[i,,] <-
                as.matrix(naive_effects)

        }

        ## quantiles of simulation distribution

        strata.quantiles = sapply(quantiles, function(q){
            apply(strata_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(strata.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

        effects.quantiles = sapply(quantiles, function(q){
            apply(effects_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(effects.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

        naive_strata.quantiles = sapply(quantiles, function(q){
            apply(naive_strata_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(naive_strata.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

        naive_effects.quantiles = sapply(quantiles, function(q){
            apply(naive_effects_MC,c(2,3),quantile,q,na.rm=T)
        },simplify="array")
        dimnames(naive_effects.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

        if (sensitivity){

            sens_strata.quantiles = sapply(quantiles, function(q){
                apply(sens_strata_MC,c(2,3),quantile,q,na.rm=T)
            },simplify="array")
            dimnames(sens_strata.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

            sens_effects.quantiles = sapply(quantiles, function(q){
                apply(sens_effects_MC,c(2,3),quantile,q,na.rm=T)
            },simplify="array")
            dimnames(sens_effects.quantiles)[[3]] = sprintf("quantile_%0.3f",quantiles)

        }

    }

    rownames(bounds$min_ca) = colnames(bounds$min_ca) = NULL
    rownames(bounds$max_ca) = colnames(bounds$max_ca) = NULL
    names(dimnames(y_sa)) = c("C", "A")
    rownames(effects.ind) = NULL

    out = list(strata = bounds$strata,
               effects = bounds$effects,
               naive = bounds$naive,
               n = nrow(data),
               n_sc_free = n_sc,
               n_sa_force = n_sa_force
               )

    if (posterior){
        out$posterior = list(
            strata.quantiles = strata.quantiles,
            effects.quantiles = effects.quantiles,
            naive_strata.quantiles = naive_strata.quantiles,
            naive_effects.quantiles = naive_effects.quantiles
        )
    }

    if (sensitivity){
        out$sensitivity = bounds$sensitivity
        if (posterior){
            out$posterior$sens_strata.quantiles = sens_strata.quantiles
            out$posterior$sens_effects.quantiles = sens_effects.quantiles
        }
    }

    return(out)
}
