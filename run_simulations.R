'%.%' <- function(x,y) paste(x,y,sep='')

## data manipulation
library(reshape2)
library(plyr)

## math
library(MASS)    # for mvrnorm()
library(gtools)  # for rdirichlet()
library(mlogit)

d.raw <- read.csv('./data/Spring2014omnibus_raw.csv', stringsAsFactors = FALSE)

source('./functions/bounds_frechet_function.R')
source('./functions/lll.R')

if (!dir.exists('results')){
  dir.create('results')
}

choice_divergences <- c(rep(0, 3), seq(0, 1, 1/3))
outcome_divergences <- c(seq(0, 1, 1/3), rep(0, 3))
draw_start <- 1
draw_stop <- 500



#####################
## preprocess data ##
#####################

## media assignment variables: convert NA/1 to 0/1
d.raw$fox[is.na(d.raw$fox)] <- 0
d.raw$msnbc[is.na(d.raw$msnbc)] <- 0
d.raw$entertainment[is.na(d.raw$entertainment)] <- 0
d.raw$med_choice[is.na(d.raw$med_choice)] <- 0  # 0 if no choice given, else 1-4

## construct 'stated preference' variable
d.raw$med_pref[which(d.raw$med_pref == 4)] <- 3
d.raw$med_pref[which(d.raw$med_pref == 0)] <- NA

## construct 'actual choice' variable
##   (pool both entertainment options in free-choice arm,
##    already pooled in forced choice arm)
d.raw$med_choice[which(d.raw$med_choice == 4)] <- 3
d.raw$med_choice[which(d.raw$med_choice == 0)] <- NA

## construct 'received treatment' assignment variable
##   (A=c if D=0, A~cat(1/3,1/3,1/3) if D=1)
##   (fox=1, msnbc=2, entertainment=3)
d.raw$med_assign <-
  ifelse(d.raw$forcedchoice == 1,
         d.raw$fox + 2*d.raw$msnbc + 3*d.raw$entertainment,
         d.raw$med_choice
         )

## construct party id variable
d.raw$pid <- NA
d.raw$pid[d.raw$party1 == 1] <- -1  # considers self a democrat
d.raw$pid[d.raw$party1 == 2] <- 1   # considers self a republican
d.raw$pid[d.raw$party4 == 1] <- 1   # neither but closer to rep
d.raw$pid[d.raw$party4 == 2] <- -1  # neither but closer to rep
d.raw$pid[d.raw$party4 == 3] <- 0   # closer to neither rep or dem

## construct pro/counter-attitudinal media variable
##   pro-attitudinal = 1      (fox for reps, or msnbc for dems)
##   counter-attitudinal = 2  (msnbc for reps, or fox for dems)
##   entertainment = 3        (for both parties)
med_labels <- c('pro-attitudinal','counter-attitudinal','entertainment')
d.raw <- d.raw[-which(d.raw$pid == 0 | is.na(d.raw$pid)),] # drop independents
d.raw$med_pref.att <-
  ifelse((d.raw$med_pref == 1 & d.raw$pid == 1) | (d.raw$med_pref == 2 & d.raw$pid == -1), 1,
         ifelse((d.raw$med_pref == 2 & d.raw$pid == 1) | (d.raw$med_pref == 1 & d.raw$pid == -1), 2, 3)
         )
d.raw$med_choice.att <-
  ifelse((d.raw$med_choice == 1 & d.raw$pid == 1) | (d.raw$med_choice == 2 & d.raw$pid == -1), 1,
         ifelse((d.raw$med_choice == 2 & d.raw$pid == 1) | (d.raw$med_choice == 1 & d.raw$pid == -1), 2, 3)
         )
d.raw$med_assign.att <-
  ifelse((d.raw$med_assign == 1 & d.raw$pid == 1) | (d.raw$med_assign == 2 & d.raw$pid == -1), 1,
         ifelse((d.raw$med_assign == 2 & d.raw$pid == 1) | (d.raw$med_assign == 1 & d.raw$pid == -1), 2, 3)
         )

## construct media sentiment index and rescale to [0,1]
##   these are eight questions on perception of media  (7-point scales)
##   name of variable is lower end of scale
##   first make direction consistent between questions
d.raw$perceive.index <- rowSums(cbind(d.raw$fair_4,
                                      d.raw$friendly_4,
                                      d.raw$good_4,
                                      d.raw$quarrel_4 * -1 + 8,
                                      d.raw$balanced_4,
                                      d.raw$oneside_4* -1 + 8,
                                      d.raw$american_4,
                                      d.raw$accurate_4
                                      ))
## then flip so larger numbers mean 'media is good'
##   for now we'll keep integer values for simulation
##   instead of rescaling like in main analysis (convert later)
d.raw$perceive.index <- -(d.raw$perceive.index - 8) + 48

## make clean data for analysis with [0,1] outcome:
##   media sentiment index
d <- data.frame(S = d.raw$med_pref.att,
                C = d.raw$med_choice.att,
                A = d.raw$med_assign.att,
                D = d.raw$forcedchoice
                )
d$Y <- d.raw$perceive.index
drop <- which(is.na(d$S) |
                is.na(d$D) |
                is.na(d$A) |
                is.na(d$Y) |
                (is.na(d$C) & d$D == 0)
              )
if (length(drop) > 0){
  d <- d[-drop,]
}
rm(drop)

## prepare to estimate effect of changing between these treatment values
J <- 3
effects.mat <- matrix(FALSE, J, J)
dimnames(effects.mat) <- list(ctrl = NULL,
                              treat = NULL)
effects.mat[2, 1] <- TRUE
effects.mat[3, 1] <- TRUE
effects.mat[3, 2] <- TRUE



########################
## summary statistics ##
########################

# sum(1{S=s,C=c})
n_sc = table(S=d$S,C=d$C)
rownames(n_sc) = "S." %.% 1:J
colnames(n_sc) = "C." %.% 1:J

# sum(1{S=s,A=a})
n_sa_pool = table(S=d$S,A=d$A)
rownames(n_sa_pool) = "S." %.% 1:J
colnames(n_sa_pool) = "A." %.% 1:J

# sum(1{S=s,A=a,D=1})
n_sa_force = table(S=d$S[d$D==1],A=d$A[d$D==1])
rownames(n_sa_force) = "S." %.% 1:J
colnames(n_sa_force) = "A." %.% 1:J

# Pr(S=s)
p_s = table(d$S) / nrow(d)

# Pr(C=c|S=s)
p_c.s = t(sapply(1:J,function(s){
  out = table(d[d$S==s & d$D==0,"C"]) / length(d[d$S==s & d$D==0,"C"])
  names(out) = "C." %.% names(out)
  return(out)
}))
rownames(p_c.s) = "S." %.% 1:J

# Pr(S=s,C=c)
p_sc_free = n_sc / sum(n_sc)

p_sc_pool = p_s %o% rep(1,J) * p_c.s
dimnames(p_sc_pool) = dimnames(p_sc_free)

p_c = colSums(p_sc_pool)

p_s.c = p_sc_pool / (rep(1,J) %o% p_c)



### observed outcomes ###

# from forced condition: E[Y|S=s,A=a,D=1]
y_sa = t(sapply(1:J,function(s){
  out = sapply(1:J,function(a){
    mean(d$Y[which(d$S==s & d$A==a & d$D==1)],na.rm=T)
  })
  names(out) = "A." %.% 1:J
  return(out)
}))
rownames(y_sa) = "S." %.% 1:J

# from free condition: E[Y|S=s,C=c,A=c,D=0]
pi_scc = t(sapply(1:J,function(s){
  out = sapply(1:J,function(c){
    mean(d$Y[which(d$S==s & d$C==c & d$D==0)],na.rm=T)
  })
  names(out) = "C." %.% 1:J
  return(out)
}))
rownames(pi_scc) = "S." %.% 1:J



##################################################
# create DGP that is consistent with pilot study #
##################################################

sc.comb = as.matrix(expand.grid(c=1:J,s=1:J)[,2:1])
sa.comb = sc.comb; colnames(sa.comb) = c("s","a")
sca.comb = as.matrix(expand.grid(a=1:J,c=1:J,s=1:J)[,3:1])

make_sim_data = function(n, choice_divergence = 0, outcome_divergence = 0){

  u = choice_divergence
  v = outcome_divergence

  ## simulation parameter u in [0,1] determines how informative the
  ##   stated choices:
  ##     u=0 means that stated choices are as informative as in real data
  ##     u=1 means that stated choices are totally uninformative and
  ##         actual choices are uniformly distributed regardless of statement
  observed_sc = p_c.s
  uninformative_sc = matrix(1/J, J, J, dimnames = dimnames(p_c.s))
  p_c.s_simpop = u * uninformative_sc + (1-u) * observed_sc

  ## fY_sca[s,c,a,] is vector of Pr( Y(a)=y | S=s, C=c )
  ##   where y = 0, ..., 48
  fY_sca_simpop = array(NA, dim=c(J,J,J,49))
  dimnames(fY_sca_simpop) = list(S = "S." %.% 1:J,
                                 C = "C." %.% 1:J,
                                 A = "A." %.% 1:J,
                                 p_Y = "Y." %.% 0:48
                                 )

  # plug in data from free-choice condition
  for (row in 1:nrow(sc.comb)){
    s = sc.comb[row,1]
    c = sc.comb[row,2]
    out = table(d$Y[d$S==s & d$C==c & d$A==c & d$D==0])
    names(out) = "Y." %.% names(out)
    out = out["Y." %.% 0:48]
    out[is.na(out)] = 0
    fY_sca_simpop[s,c,c,] = out / sum(out)
  }

  ## fY_sa[s,a,] is vector of Pr( Y(a)=y | S=s )
  ##   where y = 0, ..., 48
  fY_sa = array(NA, dim=c(J,J,49))
  dimnames(fY_sa) = list(S = "S." %.% 1:J,
                         A = "A." %.% 1:J,
                         p_Y = "Y." %.% 0:48
                         )
  # plug in data from forced-choice condition
  for (row in 1:nrow(sa.comb)){
    s = sa.comb[row,1]
    a = sa.comb[row,2]
    out = table(d$Y[d$S==s & d$A==a & d$D==1])
    names(out) = "Y." %.% names(out)
    out = out["Y." %.% 0:48]
    out[is.na(out)] = 0
    fY_sa[s,a,] = out / sum(out)
  }

  ## Observed forced-choice f( Y(a) | S=s ) is mixture of 3 components,
  ##   f( Y(a) | S=s, C=c),
  ##   where f( Y(a) | S=s, C=a ) is copied from free choice condition data.
  ## Simulation parameter outcome_divergence (v in [0, 1]) determines
  ##   how divergent the remaining two components are:
  ##     v=0 means that they are identical in distribution and
  ##     v=1 means that max of one <= min of other

  for (row in 1:nrow(sa.comb)){
    s = sa.comb[row,1]
    a = sa.comb[row,2]

    forced_mixture = fY_sa[s,a,]
    free_component = fY_sca_simpop[s,a,a,] * p_c.s_simpop[s,a]

    ## fix issues arising from the fact that our data is a finite sample
    unallocated = pmax(forced_mixture - free_component, 0)
    unallocated = unallocated / sum(unallocated) * (1 - p_c.s_simpop[s,a])

    c1 = c(2, 3, 1)[c(2, 3, 1) != a][1]
    c2 = (1:J)[!1:J %in% c(a, c1)]
    ## most balanced allocation:
    ##   divide unexplained distr evenly among choice groups by their size
    most.balanced.allocation =
      unallocated * (p_c.s_simpop[s,c1] / (1 - p_c.s_simpop[s,a]))
    ## most extreme allocation:
    ##   place portion of distr with lowest values into first choice group
    ##   (up to its size) and highest into other choice group
    unallocated.cumulative = cumsum(unallocated)
    most.extreme.allocation = unallocated
    most.extreme.allocation[unallocated.cumulative > p_c.s_simpop[s,c1]] = 0
    split.mass.at.this.y = which(unallocated.cumulative > p_c.s_simpop[s,c1])[1]
    most.extreme.allocation[split.mass.at.this.y] =
      p_c.s_simpop[s,c1] - c(0, unallocated.cumulative)[split.mass.at.this.y]

    ## select allocation by smoothing between extreme and balanced
    allocation = v * most.extreme.allocation + (1-v) * most.balanced.allocation
    fY_sca_simpop[s,c1,a,] = allocation / p_c.s_simpop[s,c1]
    fY_sca_simpop[s,c2,a,] = (unallocated - allocation) / p_c.s_simpop[s,c2]

  }

  ## stated & true prefs
  p_sc_simpop = p_s %o% rep(1,J) * p_c.s_simpop
  p_c_simpop = colSums(p_sc_simpop)
  p_s.c_simpop = p_sc_simpop / (rep(1,J) %o% p_c_simpop)

  ## Pr(Y(a)<=y | S=s)
  F_Y_sa = t(sapply(1:J,function(s){
    sapply(1:J,function(a){
      stepfun(
        x = 0:48,
        y = c(cumsum(p_c.s_simpop[s,] %*% fY_sca_simpop[s,,a,]),
              1
              )
      )
    })
  }))
  dimnames(F_Y_sa) = list(S=NULL,A=NULL)

  ## Pr(Y(c)<=y | S=s,C=c)
  F_Y_scc = t(sapply(1:J,function(s){
    sapply(1:J,function(c){
      stepfun(
        x = 0:48,
        y = c(cumsum(fY_sca_simpop[s,c,c,]), 1)
      )
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
        out = (F_Y_sa[[s,a]](y) - F_Y_scc[[s,a]](y) * p_c.s_simpop[s,a]) / (1 - p_c.s_simpop[s,a])
        return(out)
      }
      environment(f_sa) = env_sa
      attr(f_sa,"knots") = sort(unique(c(0,48,knots(F_Y_scc[[s,a]]),knots(F_Y_sa[[s,a]]))))
      attr(f_sa,"values") = sapply(attr(f_sa,"knots"),f_sa)
      return(f_sa)
    })
  }))
  dimnames(F_diff_Y_sa) = list(S=NULL,A=NULL)

  ## sampling a new dataset

  n_sc_true = matrix(rmultinom(1,n,as.numeric(p_sc_simpop)),J,J,
                     dimnames=list(S="S." %.% 1:J,C="C." %.% 1:J))

  data = do.call(rbind,sapply(1:J,function(c){
    do.call(rbind,sapply(1:J,function(s){
      data.frame(S=rep(s,n_sc_true[s,c]),
                 C=rep(c,n_sc_true[s,c]))
    },simplify=FALSE))
  },simplify=FALSE))

  data$D = sample(rep(0:1,each=n/2))

  data$A[data$D==0] = data$C[data$D==0]
  data$A[data$D==1] = sample(1:J,n/2,replace=T)

  data$Y = NA
  for (s in 1:J){
    for (c in 1:J){
      for (a in 1:J){
        data$Y[data$S==s & data$C==c & data$A==a] =
          sample(0:48, size=sum(data$S==s & data$C==c & data$A==a),
                 prob=fY_sca_simpop[s,c,a,],
                 replace=T)
      }
    }
  }

  data$C[data$D==1] = NA

  return(list(data=data,
              fY_sca_simpop = fY_sca_simpop,
              p_sc_simpop = p_sc_simpop,
              p_c_simpop = p_c_simpop,
              p_s.c_simpop = p_s.c_simpop,
              p_c.s_simpop = p_c.s_simpop,
              F_Y_sa = F_Y_sa,
              F_Y_scc = F_Y_scc,
              F_diff_Y_sa = F_diff_Y_sa
              ))

}



###############################################
# draw simulated dataset and calculate bounds #
###############################################

n_dir_MC <- 100
n_norm_MC <- 10

for (i in seq_along(choice_divergences)){

  choice_divergence <- choice_divergences[i]
  outcome_divergence <- outcome_divergences[i]

  pop.path <- file.path('results',
                        sprintf('choice%0.2f_outcome%0.2f_pop.rds',
                                choice_divergence,
                                outcome_divergence
                                )
                        )

  if (!exists(pop.path)){

    set.seed(02139)
    pop <- make_sim_data(1000,
                         choice_divergence = choice_divergence,
                         outcome_divergence = outcome_divergence
                         )
    pop$data <- NULL

    ## E[ Y(a) | S=s, C=c ]
    pop$EY_sca <-
      apply(pop$fY_sca_simpop, c(1,2,3), function(x) sum(x * 0:48))
    ## E[ Y(a) | C=c ]
    pop$EY_ca <-
      apply(pop$EY_sca * (pop$p_s.c_simpop %o% rep(1, J)), c(2, 3), sum)

    pop$bounds <- bounds.frechet.calc(
      pop$F_Y_sa,
      pop$F_Y_scc,
      pop$F_diff_Y_sa,
      pop$p_s.c_simpop,
      pop$p_c.s_simpop,
      n_sa_force = NA,
      n_sc= NA,
      y.min = 0,
      y.max = 48,
      posterior = FALSE,
      sensitivity = FALSE
    )

    pop$sens <- sens.calc(
      pop$bounds$min_ca,
      pop$bounds$max_ca,
      pop$bounds$naive_sa,
      pop$p_c.s_simpop,
      pop$p_s.c_simpop,
      effects.mat,
      seq(0, 12, 1),
      verbose = 0
    )

    saveRDS(pop, pop.path)

  }
}



#########################
## bootstrap functions ##
#########################

bounds.frechet.boot <- function(){

  ## bootstrapping
  bounds.boot.list = lapply(1:500, function(j){
    cat('bounds boot', j, '\n')
    tryCatch({
      bounds.frechet(d.sim$data[sample.int(nrow(d.sim$data), replace = TRUE),],
                     effects.mat=effects.mat,
                     posterior=FALSE,
                     sensitivity = TRUE,
                     rhos = seq(0, 12, 1)
                     )
    },
    error = function(e){
      NA
    })
  })

  ## process bootstrap results and convert to ci

  failed_bootstrap_draw = sapply(bounds.boot.list, length) == 1
  bounds.boot.list = bounds.boot.list[!failed_bootstrap_draw]

  ## naive strata
  naive_strata.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x$naive[['strata']])
    },
    simplify = 'array'
  )
  naive_strata.boot <- as.data.frame(apply(naive_strata.boot.arr, 1:2, quantile, .025))
  colnames(naive_strata.boot)[match('naive', colnames(naive_strata.boot))] <-
    'min_cilo'
  naive_strata.boot$max_cihi <- apply(naive_strata.boot.arr[,'naive',], 1, quantile, .975)

  ## naive effects
  naive_effects.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x$naive[['effects']])
    },
    simplify = 'array'
  )
  naive_effects.boot <- as.data.frame(apply(naive_effects.boot.arr, 1:2, quantile, .025))
  colnames(naive_effects.boot)[match('naive', colnames(naive_effects.boot))] <-
    'min_cilo'
  naive_effects.boot$max_cihi <- apply(naive_effects.boot.arr[,'naive',], 1, quantile, .975)

  ## bounds on strata
  strata.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x[['strata']])
    },
    simplify = 'array'
  )
  strata.boot <- as.data.frame(apply(strata.boot.arr, 1:2, quantile, .025))
  strata.boot$max <- apply(strata.boot.arr[,'max',], 1, quantile, .975)
  colnames(strata.boot)[match(c('min', 'max'), colnames(strata.boot))] <-
    c('min_cilo', 'max_cihi')

  ## bounds on effects
  effects.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x[['effects']])
    },
    simplify = 'array'
  )
  effects.boot <- as.data.frame(apply(effects.boot.arr, 1:2, quantile, .025))
  effects.boot$max <- apply(effects.boot.arr[,'max',], 1, quantile, .975)
  colnames(effects.boot)[match(c('min', 'max'), colnames(effects.boot))] <-
    c('min_cilo', 'max_cihi')

  ## sensitivity on strata
  sens_strata.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x$sensitivity[['strata']])
    },
    simplify = 'array'
  )
  sens_strata.boot <-
    as.data.frame(apply(sens_strata.boot.arr, 1:2, quantile, .025, na.rm = TRUE))
  sens_strata.boot$max <-
    apply(sens_strata.boot.arr[,'max',], 1, quantile, .975, na.rm = TRUE)
  colnames(sens_strata.boot)[match(c('min', 'max'), colnames(sens_strata.boot))] <-
    c('min_cilo', 'max_cihi')
  sens_strata.boot$failed_lp <-
    apply(sens_strata.boot.arr[,'min',], 1, function(x) mean(is.na(x)))

  ## sensitivity on effects
  sens_effects.boot.arr = sapply(
    bounds.boot.list,
    function(x){
      as.matrix(x$sensitivity[['effects']])
    },
    simplify = 'array'
  )
  sens_effects.boot <-
    as.data.frame(apply(sens_effects.boot.arr, 1:2, quantile, .025, na.rm = TRUE))
  sens_effects.boot$max <-
    apply(sens_effects.boot.arr[,'max',], 1, quantile, .975, na.rm = TRUE)
  colnames(sens_effects.boot)[match(c('min', 'max'), colnames(sens_effects.boot))] <-
    c('min_cilo', 'max_cihi')
  sens_effects.boot$failed_lp <-
    apply(sens_effects.boot.arr[,'min',], 1, function(x) mean(is.na(x)))

  return(list(
    naive_strata.ci = naive_strata.boot,
    naive_effects.ci = naive_effects.boot,
    strata.ci = strata.boot,
    effects.ci = effects.boot,
    sens_strata.ci = sens_strata.boot,
    sens_effects.ci = sens_effects.boot,
    failed_boot = mean(failed_bootstrap_draw)
  ))

}



lll.boot <- function(){

  ## don't run LLL estimator for large simulation sample sizes
  ##   (takes too long and performs poorly anyways)
  if (nrow(d.sim$data) > 3000){
    return(NULL)
  }

  lll.est <- lll(Y = d.sim$data$Y,
                 D = d.sim$data$D,
                 A = d.sim$data$A,
                 X1 = factor(d.sim$data$S),
                 X2 = factor(d.sim$data$S),
                 y.model = "normal",
                 tol = 0.001
                 )$strata

  lll.strata <- bounds.est$strata[,c('C', 'A')]
  lll.strata$est <- lll.est[as.matrix(lll.strata[,c('C', 'A')])]

  lll.effects <- bounds.est$effects[,c('C', 'treat', 'ctrl')]
  lll.effects$est <-
    lll.est[as.matrix(lll.effects[,c('C', 'treat')])] -
    lll.est[as.matrix(lll.effects[,c('C', 'ctrl')])]

  ## bootstrapping
  lll.boot.list = lapply(1:500, function(j){
    cat('lll boot', j, '\n')
    tryCatch({
      d.sim.boot <- d.sim$data[sample.int(nrow(d.sim$data), replace = TRUE),]
      lll(Y = d.sim.boot$Y,
          D = d.sim.boot$D,
          A = d.sim.boot$A,
          X1 = factor(d.sim.boot$S),
          X2 = factor(d.sim.boot$S),
          y.model = "normal",
          tol = 0.001
          )$strata
    },
    error = function(e){
      NA
    })
  })

  ## process bootstrap results and convert to ci
  failed_bootstrap_draw <- sapply(lll.boot.list, length) == 1

  lll.boot.strata <- sapply(
    lll.boot.list[!failed_bootstrap_draw],
    function(lll.boot){
      lll.boot[as.matrix(lll.strata[,c('C', 'A')])]
    })
  lll.strata$ci_lo <- apply(lll.boot.strata, 1, quantile, .025)
  lll.strata$ci_hi <- apply(lll.boot.strata, 1, quantile, .975)

  lll.boot.effects <- sapply(
    lll.boot.list[!failed_bootstrap_draw],
    function(lll.boot){
      lll.boot[as.matrix(lll.effects[,c('C', 'treat')])] -
        lll.boot[as.matrix(lll.effects[,c('C', 'ctrl')])]
    })
  lll.effects$ci_lo <- apply(lll.boot.effects, 1, quantile, .025)
  lll.effects$ci_hi <- apply(lll.boot.effects, 1, quantile, .975)

  return(list(strata = lll.strata,
              effects = lll.effects,
              boot.fail = mean(failed_bootstrap_draw)
              ))

}



#####################
## run simulations ##
#####################

for (param in seq_along(choice_divergences)){

  choice_divergence <- choice_divergences[param]
  outcome_divergence <- outcome_divergences[param]

  for (i in draw_start:draw_stop){
    for (n in c(500, 1000, 3000, 10000, 50000)){

      output.path <- file.path(
        'results',
        sprintf('choice%0.2f_outcome%0.2f_n%s_draw%s_ci.rds',
                choice_divergence,
                outcome_divergence,
                n,
                i
                )
      )
      if (file.exists(output.path)){
        next
      }

      set.seed(i)
      bounds.sim <- tryCatch({

        d.sim <- make_sim_data(n, choice_divergence, outcome_divergence)

        bounds <- bounds.frechet(d.sim$data,
                                 effects.mat=effects.mat,
                                 posterior=TRUE,
                                 n_dir_MC = n_dir_MC,
                                 n_norm_MC = n_norm_MC,
                                 hpd = TRUE,
                                 alpha = .95,
                                 sensitivity = TRUE,
                                 rhos = seq(0, 12, 1)
                                 )

        bounds.est <- bounds[-match('posterior', names(bounds))]
        bounds.post <- bounds$posterior
        bounds.boot <- bounds.frechet.boot()
        lll.est <- lll.boot()

        list(bounds.est = bounds.est,
             bounds.post = bounds.post,
             bounds.boot = bounds.boot,
             lll.est = lll.est
             )

      },
      error = function(e){
        browser()
        e
      })

      saveRDS(bounds.sim, output.path)

    }
  }
}

