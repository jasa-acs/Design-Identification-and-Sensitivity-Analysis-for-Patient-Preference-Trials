'%.%' <- function(x,y) paste(x,y,sep='')

## data manipulation
library(reshape2)
library(plyr)

## math
library(binom)

## output
library(ggplot2)
library(xtable)


sim.dir <- './results'

plots.dir <- './plots'
dir.create(plots.dir)

tables.dir <- './tables'
dir.create(tables.dir)



###########
## style ##
###########

red_dark <- '#A53F4F'
red_medium <- '#A31F34'
red_light <- '#A9606C'

blue_dark <- '#0B356F'
blue_medium <- '#315485'
blue_lightest <- '#8F9AAA'

grey_light<- '#C2C0BF'
grey_dark <- '#8A8B8C'



#################################
## load results of simulations ##
#################################

## treatment effects to compute
J <- 3
effects.mat <- matrix(FALSE, J, J)
dimnames(effects.mat) <- list(ctrl = NULL,
                              treat = NULL)
effects.mat[2, 1] <- TRUE
effects.mat[3, 1] <- TRUE
effects.mat[3, 2] <- TRUE

## list results and parse simulation parameters from filename
results.paths <- list.files(sim.dir, full.names = TRUE)
choice.divergences <- c(rep(0, 3), seq(0, 1, 1/3))
outcome.divergences <- c(seq(0, 1, 1/3), rep(0, 3))
sample.sizes <-
  as.numeric(gsub('.*n(\\d+)_.*', '\\1', basename(results.paths)))

## load results
effects.pool <-
  ldply(seq_along(choice.divergences), function(param){

    choice.divergence <- choice.divergences[param]
    outcome.divergence <- outcome.divergences[param]

    ## population distribution
    truth <- readRDS(sprintf(
      'results/choice%0.2f_outcome%0.2f_pop.rds',
      choice.divergence,
      outcome.divergence
    ))
    effects.truth <- cbind(
      n = Inf,
      type = 'truth',
      truth$bounds$effects
    )
    effects.truth$est <-
      truth$EY_ca[as.matrix(effects.truth[,c('C', 'treat')])] -
      truth$EY_ca[as.matrix(effects.truth[,c('C', 'ctrl')])]

    sens.truth <- cbind(
      ## n = Inf,
      ## type = 'truth',
      truth$sens$effects
    )
    sens.truth$true <-
      truth$EY_ca[as.matrix(sens.truth[,c('C', 'treat')])] -
      truth$EY_ca[as.matrix(sens.truth[,c('C', 'ctrl')])]

    ## sampling distribution of bounds estimates & ci
    ##   at a particular param setting
    effects.sampdistr <- ldply(
      sort(unique(na.omit(sample.sizes))),
      function(sample.size){
        out <- ldply(
          results.paths[grep(
            sprintf('choice%0.2f_outcome%0.2f_n%s',
                    choice.divergence,
                    outcome.divergence,
                    sample.size
                    ),
            results.paths,
            fixed = TRUE
          )],
          function(result.path){
            result <- tryCatch(readRDS(result.path),
                               error = function(e){
                                 e
                               })
            if ('error' %in% class(result)){
              return(data.frame())
            } else {
              effects.post <- cbind(
                type = 'posterior',
                result$bounds.est$effects,
                result$bounds.post$effects.ci[,c('min_cilo', 'max_cihi')]
              )
              colnames(effects.post) <- gsub('(min|max)_', '', colnames(effects.post))
              effects.boot <- cbind(
                type = 'bootstrap',
                result$bounds.est$effects,
                result$bounds.boot$effects.ci[,c('min_cilo', 'max_cihi')],
                fail = result$bounds.boot$failed_boot
              )
              colnames(effects.boot) <- gsub('(min|max)_', '', colnames(effects.boot))
              effects.naive.post <- cbind(
                type = 'naive',
                result$bounds.est$naive$effects,
                result$bounds.post$naive_effects.ci[,c('min_cilo', 'max_cihi')]
              )
              colnames(effects.naive.post) <- gsub('naive', 'est', colnames(effects.naive.post))
              colnames(effects.naive.post) <- gsub('(min|max)_', '', colnames(effects.naive.post))
              ## effects.naive.boot <- cbind(
              ##     type = 'naive.bootstrap',
              ##     result$bounds.est$naive$effects,
              ##     result$bounds.boot$naive_effects.ci[,c('min_cilo', 'max_cihi')]
              ## )
              ## colnames(effects.naive.boot) <- gsub('naive', 'est', colnames(effects.naive.boot))
              ## colnames(effects.naive.boot) <- gsub('(min|max)_', '', colnames(effects.naive.boot))
              if (!is.null(result$lll.est)){
                effects.lll <- cbind(
                  type = 'lll',
                  result$lll.est$effects
                )
                colnames(effects.lll) <- gsub('_', '', colnames(effects.lll))
              }
              out <- cbind(
                draw = as.numeric(gsub('.*draw(\\d+)_.*', '\\1', result.path)),
                rbind.fill(effects.post,
                           effects.boot,
                           effects.naive.post,
                           ## effects.naive.boot,
                           if (!is.null(result$lll.est)){
                             effects.lll
                           }
                           )
              )
              return(out)
            }
          })
        out <- cbind(n = sample.size, out)
        return(out)
      })

    effects.sampdistr$min.true <-
      effects.sampdistr$max.true <-
        effects.sampdistr$acte <- NA_real_
    for (i in 1:nrow(effects.truth)){
      ind <- which(
        effects.sampdistr$C == effects.truth$C[i] &
          effects.sampdistr$treat == effects.truth$treat[i] &
            effects.sampdistr$ctrl == effects.truth$ctrl[i]
      )
      effects.sampdistr$min.true[ind] <- effects.truth$min[i]
      effects.sampdistr$max.true[ind] <- effects.truth$max[i]
      effects.sampdistr$acte[ind] <- effects.truth$est[i]
    }

    effects <- rbind.fill(effects.truth, effects.sampdistr)
    effects <- cbind(choice.divergence, outcome.divergence, effects)

})

effects.pool$effect <- sprintf('treat: %s\nctrl: %s',
                               effects.pool$treat,
                               effects.pool$ctrl
                               )

effects.pool$cover.bounds <-
  ifelse(
    effects.pool$type %in% c('posterior', 'bootstrap'),
    effects.pool$cihi > effects.pool$max.true &
      effects.pool$cilo < effects.pool$min.true,
    NA
  )
effects.pool$cover.acte <-
  effects.pool$cihi > effects.pool$acte &
  effects.pool$cilo < effects.pool$acte

effects.pool$type <- as.character(effects.pool$type)

## rescale outcome to [0, 1] as in main analysis
effects.pool$min <- effects.pool$min / 48
effects.pool$max <- effects.pool$max / 48
effects.pool$est <- effects.pool$est / 48
effects.pool$min.true <- effects.pool$min.true / 48
effects.pool$max.true <- effects.pool$max.true / 48
effects.pool$acte <- effects.pool$acte / 48



##################
## analyze bias ##
##################

## quantiles <- c(.025, .1, .9, .975)
## names(quantiles) <- sprintf('quantile_%0.3f', quantiles)

C.of.interest <- 1
treat.of.interest <- 1
ctrl.of.interest <- 3

effects.pool <- effects.pool[
  effects.pool$C == C.of.interest &
    effects.pool$treat == treat.of.interest &
    effects.pool$ctrl == ctrl.of.interest,
]

bias.acte <- ddply(
  effects.pool[effects.pool$type != 'truth' &
                 grepl('lll|naive', effects.pool$type),],
  c('type', 'n', 'choice.divergence', 'outcome.divergence', 'treat', 'ctrl', 'C'),
  function(d){
    if (length(unique(d$acte)) > 1){
      stop('true acte should be same for all simulation runs')
    }
    out <- c(acte.mean = mean(d$est),
             acte.true = d$acte[1],
             bias = mean(d$est) - d$acte[1],
             se = sd(d$est) / sqrt(nrow(d))
             )
    ## out <- c(mean(d$est), quantile(d$est, quantiles))
    ## out <- out - d$acte[1]
    ## names(out) <- c('mean', names(quantiles))
    out <- as.data.frame(t(out))
    return(out)
  })
colnames(bias.acte) <- gsub('divergence', 'div', colnames(bias.acte))

bias.bounds <- ddply(
  effects.pool[effects.pool$type != 'truth' &
                 effects.pool$type == 'posterior',],
  c('type', 'n', 'choice.divergence', 'outcome.divergence', 'treat', 'ctrl', 'C'),
  function(d){
    if (length(unique(d$max.true)) > 1 || length(unique(d$min.true)) > 1){
      stop('true bounds should be same for all simulation runs')
    }
    out <- c(min.mean = mean(d$min),
             min.true = d$min.true[1],
             min.bias = mean(d$min) - d$min.true[1],
             min.se = sd(d$min) / sqrt(nrow(d)),
             max.mean = mean(d$max),
             max.true = d$max.true[1],
             max.bias = mean(d$max) - d$max.true[1],
             max.se = sd(d$max) / sqrt(nrow(d))
             )
    out <- as.data.frame(t(out))
    return(out)
  })
bias.bounds$type <- 'bounds'
colnames(bias.bounds) <- gsub('divergence', 'div', colnames(bias.bounds))

## fixed choice/outcome divergence, increasing n (sanity check):
##   mean of sampling distribution doesn't move, sampling distribution tightens
bias.acte[bias.acte$C == C.of.interest &
            bias.acte$treat == treat.of.interest &
            bias.acte$ctrl == ctrl.of.interest &
            bias.acte$choice.div == 0 &
            bias.acte$outcome.div == 0,
          ]
bias.bounds[bias.bounds$C == C.of.interest &
              bias.bounds$treat == treat.of.interest &
              bias.bounds$ctrl == ctrl.of.interest &
              bias.bounds$choice.div == 0 &
              bias.bounds$outcome.div == 0,
            ]

## ## fixed n and choice divergence, increasing outcome divergence:
## ##   bias increases for lll/naive, not for bounds
## bias.acte[bias.acte$C == C.of.interest &
##             bias.acte$treat == treat.of.interest &
##             bias.acte$ctrl == ctrl.of.interest &
##             bias.acte$n == 3000 &
##             bias.acte$choice.div == 0,
##           ]
## bias.bounds[bias.bounds$C == C.of.interest &
##               bias.bounds$treat == treat.of.interest &
##               bias.bounds$ctrl == ctrl.of.interest &
##               bias.bounds$n == 3000 &
##               bias.bounds$choice.div == 0,
##             ]

## ## fixed n and outcome divergence, increasing choice divergence:
## ##   bias increases for lll/naive, not for bounds
## bias.acte[bias.acte$C == C.of.interest &
##             bias.acte$treat == treat.of.interest &
##             bias.acte$ctrl == ctrl.of.interest &
##             bias.acte$n == 3000 &
##             bias.acte$outcome.div == 0,
##           ]
## bias.bounds[bias.bounds$C == C.of.interest &
##               bias.bounds$treat == treat.of.interest &
##               bias.bounds$ctrl == ctrl.of.interest &
##               bias.bounds$n == 3000 &
##               bias.bounds$outcome.div == 0,
##             ]


## REPLICATE: table 2

tab2.naive.cd <- bias.acte[bias.acte$C == C.of.interest &
                             bias.acte$treat == treat.of.interest &
                             bias.acte$ctrl == ctrl.of.interest &
                             bias.acte$n == 3000 &
                             bias.acte$type == 'naive' &
                             bias.acte$outcome.div == 0,
                           c('type', 'choice.div', 'bias')
                           ]
tab2.naive.cd$choice.div <- sprintf('CD=%0.2f', tab2.naive.cd$choice.div)

tab2.naive.od <- bias.acte[bias.acte$C == C.of.interest &
                             bias.acte$treat == treat.of.interest &
                             bias.acte$ctrl == ctrl.of.interest &
                             bias.acte$n == 3000 &
                             bias.acte$type == 'naive' &
                             bias.acte$choice.div == 0,
                           c('type', 'outcome.div', 'bias')
                           ]
tab2.naive.od$outcome.div <- sprintf('OD=%0.2f', tab2.naive.od$outcome.div)

tab2.bounds.cd <- bias.bounds[bias.bounds$C == C.of.interest &
                                bias.bounds$treat == treat.of.interest &
                                bias.bounds$ctrl == ctrl.of.interest &
                                bias.bounds$n == 3000 &
                                bias.bounds$outcome.div == 0,
                              c('choice.div', 'min.bias', 'max.bias')
                              ]
tab2.bounds.cd$choice.div <- sprintf('CD=%0.2f', tab2.bounds.cd$choice.div)

tab2.bounds.od <- bias.bounds[bias.bounds$C == C.of.interest &
                                bias.bounds$treat == treat.of.interest &
                                bias.bounds$ctrl == ctrl.of.interest &
                                bias.bounds$n == 3000 &
                                bias.bounds$choice.div == 0,
                              c('outcome.div', 'min.bias', 'max.bias')
                              ]
tab2.bounds.od$outcome.div <- sprintf('OD=%0.2f', tab2.bounds.od$outcome.div)

xtable(
  rbind(
    acast(
      tab2.naive.cd,
      type ~ choice.div,
      value.var = 'bias',
      ),
    acast(
      melt(tab2.bounds.od, id.vars = 'outcome.div'),
      variable ~ outcome.div
    )
  ),
  digits = 3
)

xtable(
  rbind(
    acast(
      tab2.naive.od,
      type ~ outcome.div,
      value.var = 'bias',
      ),
    acast(
      melt(tab2.bounds.cd, id.vars = 'choice.div'),
      variable ~ choice.div
    )
  ),
  digits = 3
)

tab2.lll.cd <- bias.acte[bias.acte$C == C.of.interest &
                           bias.acte$treat == treat.of.interest &
                           bias.acte$ctrl == ctrl.of.interest &
                           bias.acte$n == 3000 &
                           bias.acte$type == 'lll' &
                           bias.acte$outcome.div == 0,
                         c('type', 'choice.div', 'bias')
                         ]
tab2.lll.cd$choice.div <- sprintf('CD=%0.2f', tab2.lll.cd$choice.div)

tab2.lll.od <- bias.acte[bias.acte$C == C.of.interest &
                           bias.acte$treat == treat.of.interest &
                           bias.acte$ctrl == ctrl.of.interest &
                           bias.acte$n == 3000 &
                           bias.acte$type == 'lll' &
                           bias.acte$choice.div == 0,
                         c('type', 'outcome.div', 'bias')
                         ]
tab2.lll.od$outcome.div <- sprintf('OD=%0.2f', tab2.lll.od$outcome.div)

xtable(
  rbind(
    acast(
      tab2.lll.cd,
      type ~ choice.div,
      value.var = 'bias',
      )
  ),
  digits = 3
)

xtable(
  rbind(
    acast(
      tab2.lll.od,
      type ~ outcome.div,
      value.var = 'bias',
      )
  ),
  digits = 3
)



######################
## analyze coverage ##
######################

cvg <- ddply(
  effects.pool[effects.pool$type != 'truth',],
  c('type', 'n', 'choice.divergence', 'outcome.divergence', 'treat', 'ctrl', 'C'),
  function(d){
    test.acte <- binom.confint(x = sum(d$cover.acte, na.rm = TRUE),
                               n = length(na.omit(d$cover.acte)),
                               method = 'exact'
                               )[,c('mean', 'lower', 'upper')]
    colnames(test.acte) <- c('cvg.acte',
                             'cvg.acte.cilo',
                             'cvg.acte.cihi'
                             )
    if (d$type[1] %in% c('bootstrap', 'posterior')){
      test.bounds <- binom.confint(x = sum(d$cover.bounds, na.rm = TRUE),
                                   n = length(na.omit(d$cover.bounds)),
                                   method = 'exact'
                                   )[,c('mean', 'lower', 'upper')]
      colnames(test.bounds) <- c('cvg.bounds',
                                 'cvg.bounds.cilo',
                                 'cvg.bounds.cihi'
                                 )
      test.bounds$fail <- mean(d$fail)
      return(cbind(test.acte, test.bounds))
    } else {
      return(test.acte)
    }
  })
colnames(cvg) <- gsub('divergence', 'div', colnames(cvg))

## print(
##   dcast(
##     cvg[cvg$C == C.of.interest &
##           cvg$treat == treat.of.interest &
##           cvg$ctrl == ctrl.of.interest &
##           ## cvg$n == 3000 &
##           ## cvg$outcome.div == 0 &
##           cvg$choice.div == 0,
##         c('type', 'n', 'cvg.acte', 'outcome.div')
##         ],
##     type + n ~ outcome.div,
##     value.var = 'cvg.acte'
##   ),
##   digits = 2,
##   row.names = FALSE
## )

cvg$choice.div <- sprintf('CD=%0.2f', cvg$choice.div)
cvg$outcome.div <- sprintf('OD=%0.2f', cvg$outcome.div)

print(
  xtable(
    dcast(
      cvg[cvg$C == C.of.interest &
            cvg$treat == treat.of.interest &
            cvg$ctrl == ctrl.of.interest &
            cvg$type == 'posterior' &
            ## cvg$n == 3000 &
            ## cvg$outcome.div == 0 &
            cvg$choice.div == 'CD=0.00',
          c('n', 'cvg.bounds', 'outcome.div')
          ],
      n ~ outcome.div,
      value.var = 'cvg.bounds'
    ),
    digits = 3
  ),
  include.rownames = FALSE
)

print(
  xtable(
    dcast(
      cvg[cvg$C == C.of.interest &
            cvg$treat == treat.of.interest &
            cvg$ctrl == ctrl.of.interest &
            cvg$type == 'posterior' &
            cvg$outcome.div == 'OD=0.00'
         ,
          c('type', 'n', 'cvg.bounds', 'choice.div')
          ],
      n ~ choice.div,
      value.var = 'cvg.bounds'
    ),
    digits = 3
  ),
  include.rownames = FALSE
)

print(
  xtable(
    dcast(
      cvg[cvg$C == C.of.interest &
            cvg$treat == treat.of.interest &
            cvg$ctrl == ctrl.of.interest &
            cvg$type == 'bootstrap' &
            cvg$choice.div == 'CD=0.00',
          c('n', 'cvg.bounds', 'outcome.div')
          ],
      n ~ outcome.div,
      value.var = 'cvg.bounds'
    ),
    digits = 3
  ),
  include.rownames = FALSE
)

print(
  xtable(
    dcast(
      cvg[cvg$C == C.of.interest &
            cvg$treat == treat.of.interest &
            cvg$ctrl == ctrl.of.interest &
            cvg$type == 'bootstrap' &
            cvg$outcome.div == 'OD=0.00'
         ,
          c('type', 'n', 'cvg.bounds', 'choice.div')
          ],
      n ~ choice.div,
      value.var = 'cvg.bounds'
    ),
    digits = 3
  ),
  include.rownames = FALSE
)

