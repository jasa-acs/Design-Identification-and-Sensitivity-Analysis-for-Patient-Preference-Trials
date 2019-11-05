###############################################################
## replication file:                                         ##
## design, identification, and sensitivity analysis for ppts ##
###############################################################

## note: S_i, C_i, A_i are zero-indexed in paper but one-indexed in code

'%.%' <- function(x,y) paste(x,y,sep='')

## data manipulation
library(reshape2)
library(plyr)

## output
library(ggplot2)
library(xtable)

## math
library(MASS)    # for mvrnorm()
library(gtools)  # for rdirichlet()
library(lpSolve)

data.dir <- './data'
software.dir <- ''

d.raw <- read.csv('./data/Spring2014omnibus_raw.csv', stringsAsFactors = FALSE)

source('./functions/bounds_frechet_function.R')
source('./functions/bounds_binary_function.R')

if (!dir.exists('results')){
  dir.create('results')
}



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
##   flip all questions so larger numbers mean 'media is good'
d.raw$perceive.index <- rowSums(cbind(d.raw$fair_4,
                                      d.raw$friendly_4,
                                      d.raw$good_4,
                                      d.raw$quarrel_4 * -1 + 8,
                                      d.raw$balanced_4,
                                      d.raw$oneside_4* -1 + 8,
                                      d.raw$american_4,
                                      d.raw$accurate_4
                                      ))
d.raw$perceive.index <- (d.raw$perceive.index - 8) / -48 + 1

## make clean data for analysis with [0,1] outcome:
##   media sentiment index
d.frechet <- data.frame(S = d.raw$med_pref.att,
                        C = d.raw$med_choice.att,
                        A = d.raw$med_assign.att,
                        D = d.raw$forcedchoice
                        )
d.frechet$Y <- d.raw$perceive.index
drop <- which(is.na(d.frechet$S) |
                is.na(d.frechet$D) |
                is.na(d.frechet$A) |
                is.na(d.frechet$Y) |
                (is.na(d.frechet$C) & d.frechet$D == 0)
              )
if (length(drop) > 0){
  d.frechet <- d.frechet[-drop,]
}
rm(drop)

## make clean data for analysis with binary outcome:
##   at least somewhat likely to discuss clip with friend
d.binary <- data.frame(S = d.raw$med_pref.att,
                       C = d.raw$med_choice.att,
                       A = d.raw$med_assign.att,
                       D = d.raw$forcedchoice
                       )
d.binary$Y <- d.raw$actions_discuss <= 3  # not likely/not sure = 0, else = 1
drop <- which(is.na(d.binary$S) |
                is.na(d.binary$D) |
                is.na(d.binary$A) |
                is.na(d.binary$Y) |
                (is.na(d.binary$C) & d.binary$D == 0)
              )
if (length(drop) > 0){
  d.binary <- d.binary[-drop,]
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



######################################################
## various summary statistics reported in section 7 ##
######################################################

d.summary <- data.frame(S = d.raw$med_pref.att,
                        C = d.raw$med_choice.att,
                        A = d.raw$med_assign.att,
                        D = d.raw$forcedchoice,
                        Y.sentiment = d.raw$perceive.index,
                        Y.discuss = d.raw$actions_discuss <= 3
                        )

## REPLICATE: distr of stated prefs, p20
round(
  prop.table(
    table(S = d.summary$S)
  ),
  2
)

## REPLICATE: consistency of stated/actual preferences, p20
round(
  diag(prop.table(
    table(S = d.summary$S,
          C = d.summary$C),
    margin = 1
  )),
  2
)

## REPLICATE: table 1, upper section (free-choice)
## xtable(
round(t(
  ddply(d.summary[d.summary$D == 0,], c('S', 'A'), function(d){
    data.frame(
      strata.prop = nrow(d) / sum(d.summary$D == 0),
      sentiment = mean(d$Y.sentiment, na.rm = TRUE),
      discuss = mean(d$Y.discuss, na.rm = TRUE)
    )
  })
), 2)
## )

## REPLICATE: table 1, lower section (forced-choice)
## xtable(
round(t(
  ddply(d.summary[d.summary$D == 1,], c('S', 'A'), function(d){
    data.frame(
      strata.prop = nrow(d) / sum(d.summary$D == 0),
      sentiment = mean(d$Y.sentiment, na.rm = TRUE),
      discuss = mean(d$Y.discuss, na.rm = TRUE)
    )
  })
), 2)
## )


## REPLICATE: outcome means/sds, p20
round(mean(d.summary$Y.sentiment, na.rm = TRUE), 2)
round(sd(d.summary$Y.sentiment, na.rm = TRUE), 2)
round(mean(d.summary$Y.discuss, na.rm = TRUE), 3)



#####################################################################
## REPLICATE: left panel of figure 2, bounds for outcomes in [0,1] ##
#####################################################################

## conduct main analysis: bounds and ci for arbitrary outcomes
##   (implements procedure in sec 4.1 using
##    code from ./functions/bounds_frechet_function.R)
## this also computes sensitivity analysis
if (file.exists('./results/acte_frechet_posterior5k.rds')){
  acte.frechet <- readRDS('./results/acte_frechet_posterior5k.rds')
} else {
  set.seed(02139)
  acte.frechet <- bounds.frechet(
    data = d.frechet,
    effects.mat = effects.mat,
    posterior = TRUE,
    n_dir_MC = 500,
    n_norm_MC = 10,
    sensitivity = TRUE,
    rhos = seq(0, .25, .01),
    hpd = TRUE,
    alpha = .95
  )
  saveRDS(acte.frechet, './results/acte_frechet_posterior5k.rds')
}



### plot frechet bounds ###

plot_effects <- acte.frechet$effects

plot_effects$min.cilo <- acte.frechet$posterior$effects.ci$min_cilo
plot_effects$max.cihi <- acte.frechet$posterior$effects.ci$max_cihi

plot_effects$naive <- acte.frechet$naive$effects$naive
plot_effects$naive.cilo <- acte.frechet$posterior$naive_effects.ci$min_cilo
plot_effects$naive.cihi <- acte.frechet$posterior$naive_effects.ci$max_cihi

plot_ate <- unique(lapply(1:nrow(plot_effects),function(i){
  out <- as.numeric(plot_effects[i,c('ctrl','treat')])
  test <- t.test(d.frechet$Y[d.frechet$D==1 & d.frechet$A==out[2]],d.frechet$Y[d.frechet$D==1 & d.frechet$A==out[1]])
  out <- c(out,-diff(test$estimate),test$conf.int)
}))
plot_ate <- as.data.frame(do.call(rbind,plot_ate))
colnames(plot_ate) <- c('ctrl','treat','ate','ate.cilo','ate.cihi')
plot_ate$C <- 0

plot_effects <- rbind.fill(plot_effects,plot_ate)

plot_effects$C <- factor(plot_effects$C,levels=0:J,labels=c('pooled\n(ATE)',med_labels))
plot_effects$ctrl <- factor(plot_effects$ctrl,levels=1:J,labels= "a' = " %.% med_labels)
plot_effects$treat <- factor(plot_effects$treat,levels=1:J,labels= 'a = ' %.% med_labels)
plot_effects$treatment <- factor(plot_effects$treat %.% '\n' %.% plot_effects$ctrl)
plot_effects$treatment <- relevel(plot_effects$treatment,levels(plot_effects$treatment)[3])
plot_effects$ctrl <- plot_effects$treat <- NULL

plot_effects.melt <- melt(plot_effects,id.vars=c('C','treatment'))
plot_effects.melt$type <- ifelse(grepl('min|max',plot_effects.melt$variable),'bounds','naive')
plot_effects.melt$var <- sapply(as.character(plot_effects.melt$variable),function(x){
  switch(x,
         min='min',
         max='max',
         min.cilo='cilo',
         min.cihi=NA,
         max.cilo=NA,
         max.cihi='cihi',
         naive='point',
         naive.cilo='cilo',
         naive.cihi='cihi',
         ate='point',
         ate.cilo='cilo',
         ate.cihi='cihi'
         )
})
plot_effects.melt$variable <- NULL
plot_effects <- dcast(na.omit(plot_effects.melt), C + treatment + type ~ var)
plot_effects$type <- factor(plot_effects$type,levels=c('naive','bounds'))
plot_effects$C <-
  factor(plot_effects$C,levels=c('pooled\n(ATE)',med_labels),
         labels=c('pooled\n(ATE)','pro-\nattitudinal','counter-\nattitudinal','enter-\ntainment'),
         ordered=TRUE)

pdf('./results/fig2_FH_bounds_naive_ATE.pdf', 5, 8)
ggplot(plot_effects, aes(x = C, colour = type, linetype = type)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = cilo, ymax = cihi),
                width = .5,
                position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = min, ymax = max),
                width = .25,
                size = 1.25,
                position = position_dodge(width = .5)) +
  geom_point(aes(y = point),
             position = position_dodge(width = .5),
             size = 1.5) +
  geom_errorbar(aes(ymin = cilo, ymax = cihi),
                data = plot_effects[plot_effects$C=='pooled\n(ATE)',],
                width = .5,
                size = 1,
                linetype = 'solid') +
  geom_point(aes(y = point),
             data = plot_effects[plot_effects$C=='pooled\n(ATE)',],
             position = position_dodge(width = .5),
             size = 2.5) +
  facet_grid(treatment ~ .) +
  scale_colour_manual(values = c(naive = blue_medium, bounds = red_medium),
                      guide = FALSE) +
  scale_linetype_manual(values = c(naive = 'dashed', bounds = 'solid'),
                        guide = FALSE) +
  ylab(expression('Expected change in sentiment index, ' *
                    tau *"(a, a' | c)")) +
  xlab(expression('Among subjects who would choose ' * (C[i]))) +
  scale_y_continuous(breaks = seq(-.5, .5, .25),
                     limits = c(-.425, .425)) +
  ggtitle('Sentiment toward Media\n') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
rm(plot_effects)



############################################################
## REPLICATE: figure 3, sensitivity for outcomes in [0,1] ##
############################################################

### plot previously computed sensitivity results for frechet method ###

C.labels <- 'would choose\n' %.% med_labels
names(C.labels) <- c('1', '2', '3')

treatment.labels <- c('13' = "a = pro-attitudinal\na' = entertainment",
                      '23' = "a = counter-attitudinal\na' = entertainment",
                      '12' = "a = pro-attitudinal\na' = counter-attitudinal")

acte.frechet$sensitivity$effects$treatment <-
  factor(treatment.labels[acte.frechet$sensitivity$effects$treat %.%
                            acte.frechet$sensitivity$effects$ctrl],
         levels = treatment.labels)

## naive estimates (should be about the same as normal posterior)
acte.frechet$naive.ci <- as.data.frame(acte.frechet$naive$effects)
acte.frechet$naive.ci$min <-
  acte.frechet$posterior$naive_effects.ci$min_cilo
acte.frechet$naive.ci$max <-
  acte.frechet$posterior$naive_effects.ci$max_cihi
acte.frechet$naive.ci$treatment <-
  factor(treatment.labels[acte.frechet$naive.ci$treat %.% acte.frechet$naive.ci$ctrl],
         levels = treatment.labels)

## bounds on treatment effects
acte.frechet$posterior$effects.ci <-
  as.data.frame(acte.frechet$posterior$effects.quantiles[,, 'quantile_0.025'])
acte.frechet$posterior$effects.ci$max <-
  acte.frechet$posterior$effects.quantiles[, 'max', 'quantile_0.975']
acte.frechet$posterior$effects.ci$treatment <-
  factor(treatment.labels[acte.frechet$posterior$effects.ci$treat %.%
                            acte.frechet$posterior$effects.ci$ctrl],
         levels = treatment.labels)

## sensitivity analysis for treatment effects
acte.frechet$sensitivity$effects.ci <- as.data.frame(
  acte.frechet$posterior$sens_effects.quantiles[,, 'quantile_0.025']
)
acte.frechet$sensitivity$effects.ci$max <-
  acte.frechet$posterior$sens_effects.quantiles[, 'max', 'quantile_0.975']
acte.frechet$sensitivity$effects.ci$treatment <-
  factor(treatment.labels[acte.frechet$sensitivity$effects.ci$treat %.%
                            acte.frechet$sensitivity$effects.ci$ctrl],
         levels = treatment.labels)

## naive estimates assume rho == 0
##   and bounds equivalent to rho == Inf
##   (for plotting purposes, place bounds at maximum rho)
acte.frechet$naive$effects$rho <- 0
acte.frechet$naive$effects$treatment <-
  factor(treatment.labels[acte.frechet$naive$effects$treat %.%
                            acte.frechet$naive$effects$ctrl],
         levels = treatment.labels)

acte.frechet$naive.ci$rho <- 0
acte.frechet$naive.ci$treatment <-
  factor(treatment.labels[acte.frechet$naive.ci$treat %.%
                            acte.frechet$naive.ci$ctrl],
         levels = treatment.labels)

acte.frechet$effects$rho <-
  max(acte.frechet$sensitivity$effects$rho[
    is.finite(acte.frechet$sensitivity$effects$rho)
  ])
acte.frechet$effects$treatment <-
  factor(treatment.labels[acte.frechet$effects$treat %.%
                            acte.frechet$effects$ctrl],
         levels = treatment.labels)

acte.frechet$posterior$effects.ci$rho <-
  max(acte.frechet$sensitivity$effects$rho[
    is.finite(acte.frechet$sensitivity$effects$rho)
  ])
acte.frechet$posterior$effects.ci$treatment <-
  factor(treatment.labels[acte.frechet$posterior$effects.ci$treat %.%
                            acte.frechet$posterior$effects.ci$ctrl],
         levels = treatment.labels)

pdf('./results/fig3_FH_sensitivity.pdf', 8, 8)
ggplot(acte.frechet$sensitivity$effects[
  is.finite(acte.frechet$sensitivity$effects$rho),], aes(x=rho)) +
  geom_ribbon(aes(ymin=min, ymax=max),
              data=acte.frechet$sensitivity$effects.ci[
                is.finite(acte.frechet$sensitivity$effects.ci$rho),],
              fill='#80808080') +
  geom_ribbon(aes(ymin=min, ymax=max),
              fill='#40404080') +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'naive'),
                data = acte.frechet$naive.ci,
                width = .02, size = .5, linetype = 'dashed') +
  geom_point(aes(y = naive, colour = 'naive'),
             data = acte.frechet$naive$effects) +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'bounds'),
                data = acte.frechet$posterior$effects.ci,
                width = .015, size = .5) +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'bounds'),
                data = acte.frechet$effects,
                width = .01, size = 1) +
  facet_grid(treatment ~ C, labeller = labeller(C = C.labels)) +
  scale_color_manual(values = c(naive = blue_medium, bounds = red_medium),
                     guide = FALSE) +
  xlab(expression('Sensitivity parameter ' * (rho))) +
  ylab(expression('Expected change in sentiment index, ' * tau("a, a'| c"))) +
  scale_x_continuous(breaks = seq(0, .25, .1)) +
  scale_y_continuous(breaks=seq(-.5,.5,.25), limits = c(-.425, .425)) +
  ggtitle('Sensitivity Analysis for Sentiment toward Media\n') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



### interpret sensitivity results ###

effects.ind = which(effects.mat, arr.ind=T)
colnames(effects.ind) = c('ctrl', 'treat')

zero.crossings <- expand.grid(C = 1:J,
                              treatment.ind = 1:nrow(effects.ind)
                              )
zero.crossings$ctrl <- effects.ind[,'ctrl'][zero.crossings$treatment.ind]
zero.crossings$treat <- effects.ind[,'treat'][zero.crossings$treatment.ind]
zero.crossings$treatment.ind <- NULL
zero.crossings$min.crosses.zero <-
  zero.crossings$max.crosses.zero <- NA
zero.crossings$min.ci.crosses.zero <-
  zero.crossings$max.ci.crosses.zero <- NA
zero.crossings$min.becomes.flat <-
  zero.crossings$max.becomes.flat <- NA
zero.crossings$min.ci.becomes.flat <-
  zero.crossings$max.ci.becomes.flat <- NA

## "estimated maximal difference [of rho] for any strata is 0.18", p24
##   estimated sensitivity interval equal estimated bounds at rho = 0.18
as.matrix(acte.frechet$sensitivity$strata[
  acte.frechet$sensitivity$strata$rho == .18, c('min', 'max')]
  ) -
  as.matrix(acte.frechet$sensitivity$strata[
    acte.frechet$sensitivity$strata$rho == Inf, c('min', 'max')]
    )
## but estimated sensitivity interval don't equal estimated bounds yet
##   for rho = 0.17
as.matrix(acte.frechet$sensitivity$strata[
  acte.frechet$sensitivity$strata$rho == .17, c('min', 'max')]
  ) -
  as.matrix(acte.frechet$sensitivity$strata[
    acte.frechet$sensitivity$strata$rho == Inf, c('min', 'max')]
    )

## "for most strata, differences above 0.1 can be ruled out", p24
##   for all but one choice-specific potential outcome,
##   estimated sensitivity interval equal estimated bounds at rho = 0.1
as.matrix(acte.frechet$sensitivity$strata[
  acte.frechet$sensitivity$strata$rho == .1, c('min', 'max')]
  ) - as.matrix(acte.frechet$sensitivity$strata[
    acte.frechet$sensitivity$strata$rho == Inf, c('min', 'max')]
    )



## for a partiular ACTE, determine the values of rho at which
##   - estimated sensitivity interval no longer includes zero
##   - credible intervals for sensitivity endpoints no longer include zero
for (i in 1:nrow(zero.crossings)){

  ## subset to choice/treatment of interest
  ind <- which(
    acte.frechet$sensitivity$effects$C == zero.crossings$C[i] &
      acte.frechet$sensitivity$effects$ctrl == zero.crossings$ctrl[i] &
        acte.frechet$sensitivity$effects$treat == zero.crossings$treat[i]
  )

  ## prepare to invert estimated sensitivity bounds
  ##   (goal is to find value of rho at which bounds cross zero; to do this,
  ##    first drop values in low-rho NA region and high-rho flat region)
  min.drop <- diff(c(0, acte.frechet$sensitivity$effects$min[ind])) == 0
  max.drop <- diff(c(0, acte.frechet$sensitivity$effects$max[ind])) == 0
  zero.crossings$min.becomes.flat[i] <-
    acte.frechet$sensitivity$effects$rho[ind][
      which(min.drop)[1] - 1
    ]
  zero.crossings$max.becomes.flat[i] <-
    acte.frechet$sensitivity$effects$rho[ind][
      which(max.drop)[1] - 1
    ]
  min.drop[is.na(min.drop)] <- TRUE
  max.drop[is.na(max.drop)] <- TRUE

  ## repeat prep for upper (lower) ci of upper (lower) sensitivity endpoint
  min.ci.drop <- diff(c(
    0,
    acte.frechet$posterior$sens_effects.quantiles[
      ind, 'min', 'quantile_0.025'
    ]
  )) == 0
  max.ci.drop <- diff(c(
    0,
    acte.frechet$posterior$sens_effects.quantiles[
      ind, 'max', 'quantile_0.975'
    ]
  )) == 0
  zero.crossings$min.ci.becomes.flat[i] <-
    acte.frechet$sensitivity$effects$rho[ind][
      which(min.ci.drop)[1] - 1
    ]
  zero.crossings$max.ci.becomes.flat[i] <-
    acte.frechet$sensitivity$effects$rho[ind][
      which(max.ci.drop)[1] - 1
    ]
  min.ci.drop[is.na(min.ci.drop)] <- TRUE
  max.ci.drop[is.na(max.ci.drop)] <- TRUE

  ## invert estimated sensitivity and ci
  if (length(ind[!min.drop]) > 1){
    zero.crossings$min.crosses.zero[i] <- approx(
      x = acte.frechet$sensitivity$effects$min[ind][!min.drop],
      y = acte.frechet$sensitivity$effects$rho[ind][!min.drop],
      xout = 0
    )$y
  }
  if (length(ind[!min.ci.drop]) > 1){
    zero.crossings$min.ci.crosses.zero[i] <- approx(
      x = acte.frechet$posterior$sens_effects.quantiles[
        ind, 'min', 'quantile_0.025'
      ][!min.ci.drop],
      y = acte.frechet$posterior$sens_effects.quantiles[
        ind, 'rho', 'quantile_0.025'
      ][!min.ci.drop],
      xout = 0
    )$y
  }
  if (length(ind[!max.drop]) > 1){
    zero.crossings$max.crosses.zero[i] <- approx(
      x = acte.frechet$sensitivity$effects$max[ind][!max.drop],
      y = acte.frechet$sensitivity$effects$rho[ind][!max.drop],
      xout = 0
    )$y
  }
  if (length(ind[!max.ci.drop]) > 1){
    zero.crossings$max.ci.crosses.zero[i] <- approx(
      x = acte.frechet$posterior$sens_effects.quantiles[
        ind, 'max', 'quantile_0.975'
      ][!max.ci.drop],
      y = acte.frechet$posterior$sens_effects.quantiles[
        ind, 'rho', 'quantile_0.975'
      ][!max.ci.drop],
      xout = 0
    )$y
  }
}

## REPLICATE: "we focus on the middle plot in the center row of figure 3...", p25 top para
acte.frechet$sensitivity$effects[           # rho < 0.02 is infeasible
  acte.frechet$sensitivity$effects$treat == 2 &
    acte.frechet$sensitivity$effects$ctrl == 3 &
    acte.frechet$sensitivity$effects$C == 2,]
zero.crossings[zero.crossings$treat == 2 &  # estimated sensitivity is flat for
                 zero.crossings$ctrl == 3 &   #   rho >= 0.18, and estimates/ci
                 zero.crossings$C == 2,]      #   cross zero at the given values



######################################################################
## REPLICATE: right panel of figure 2, bounds for outcomes in {0,1} ##
######################################################################

## conduct secondary analysis: bounds and ci for binary outcomes
##   (implements procedure in sec 4.2 using
##    code from ./functions/bounds_binary_function.R)
## this also computes sensitivity analysis

if (file.exists('./results/acte_binary_posterior5k.rds')){
  acte.binary <- readRDS('./results/acte_binary_posterior5k.rds')
} else {
  set.seed(02139)
  acte.binary <- bounds.binary(
    data = d.binary,
    effects.mat = effects.mat,
    posterior = TRUE,
    n_MC = 5000,
    quantiles = c(.025,.975),
    sensitivity = TRUE,
    rhos = seq(0, 0.25, .01),
    hpd = TRUE,
    alpha = .95,
    verbose = 2
  )
  saveRDS(acte.binary, './results/acte_binary_posterior5k.rds')
}



### plot binary bounds ###

med_labels <- c('pro-attitudinal','counter-attitudinal','entertainment')

plot_effects <- acte.binary$effects

plot_effects$min.cilo <- acte.binary$posterior$effects.ci$min_cilo
## plot_effects$min.cihi <- effects.cihi$min
plot_effects$max.cihi <- acte.binary$posterior$effects.ci$max_cihi
## plot_effects$max.cihi <- effects.cihi$max

plot_effects$naive <- acte.binary$naive$effects$naive
plot_effects$naive.cilo <-
  acte.binary$posterior$naive_effects.ci$min_cilo
plot_effects$naive.cihi <-
  acte.binary$posterior$naive_effects.ci$max_cihi

plot_ate <- unique(lapply(1:nrow(plot_effects),function(i){
  out <- as.numeric(plot_effects[i,c('ctrl','treat')])
  test <- t.test(d.binary$Y[d.binary$D==1 & d.binary$A==out[2]],
                 d.binary$Y[d.binary$D==1 & d.binary$A==out[1]])
  out <- c(out,-diff(test$estimate),test$conf.int)
}))

plot_ate <- as.data.frame(do.call(rbind,plot_ate))
colnames(plot_ate) <- c('ctrl','treat','ate','ate.cilo','ate.cihi')
plot_ate$C <- 0

plot_effects <- rbind.fill(plot_effects,plot_ate)

plot_effects$C <-
  factor(plot_effects$C,levels=0:J,labels=c('pooled\n(ATE)',med_labels))
plot_effects$ctrl <-
  factor(plot_effects$ctrl,levels=1:J,labels= "a' = " %.% med_labels)
plot_effects$treat <-
  factor(plot_effects$treat,levels=1:J,labels= 'a = ' %.% med_labels)
plot_effects$treatment <-
  factor(plot_effects$treat %.% '\n' %.% plot_effects$ctrl)
plot_effects$treatment <-
  relevel(plot_effects$treatment,levels(plot_effects$treatment)[3])
plot_effects$ctrl <- plot_effects$treat <- NULL

plot_effects.melt <- melt(plot_effects,id.vars=c('C','treatment'))
plot_effects.melt$type <- ifelse(grepl('min|max', plot_effects.melt$variable),
                                 'bounds',
                                 'naive')
plot_effects.melt$var <- sapply(as.character(plot_effects.melt$variable),function(x){
  switch(x,
         min='min',
         max='max',
         min.cilo='cilo',
         min.cihi=NA,
         max.cilo=NA,
         max.cihi='cihi',
         naive='point',
         naive.cilo='cilo',
         naive.cihi='cihi',
         ate='point',
         ate.cilo='cilo',
         ate.cihi='cihi'
         )
})
plot_effects.melt$variable <- NULL
plot_effects <- dcast(na.omit(plot_effects.melt), C + treatment + type ~ var)
plot_effects$type <- factor(plot_effects$type,levels=c('naive','bounds'))
plot_effects$C <-
  factor(plot_effects$C,levels=c('pooled\n(ATE)',med_labels),
         labels=c('pooled\n(ATE)',
                  'pro-\nattitudinal',
                  'counter-\nattitudinal',
                  'enter-\ntainment'),
         ordered = TRUE)

pdf('./results/fig2_binary_bounds_naive_ATE.pdf', 5, 8)
ggplot(plot_effects, aes(x = C, colour = type, linetype = type)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = cilo, ymax = cihi),
                width = .5,
                position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = min, ymax = max),
                width = .25,
                size = 1.25,
                position = position_dodge(width = .5)) +
  geom_point(aes(y = point),
             position = position_dodge(width = .5),
             size = 1.5) +
  geom_errorbar(aes(ymin = cilo, ymax = cihi),
                data = plot_effects[plot_effects$C == 'pooled\n(ATE)',],
                width = .5,
                size = 1,
                linetype = 'solid') +
  geom_point(aes(y = point),
             data = plot_effects[plot_effects$C=='pooled\n(ATE)',],
             position = position_dodge(width = .5),
             size = 2.5) +
  facet_grid(treatment ~ .) +
  scale_colour_manual(values = c(naive = blue_medium, bounds = red_medium),
                      guide = FALSE) +
  scale_linetype_manual(values = c(naive = 'dashed', bounds = 'solid'),
                        guide = FALSE) +
  ylab(expression('Expected change in probability of discussing, ' * tau *"(a, a' | c)")) +
  xlab(expression('Among subjects who would choose ' * (C[i]))) +
  scale_y_continuous(breaks = seq(-.5, .5, .25),
                     limits = c(-.65, .65)) +
  ggtitle('Discuss Story with Friends (binary)\n') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



##############################################################
## REPLICATE: figure A.1, sensitivity for outcomes in [0,1] ##
##############################################################

### plot previously computed sensitivity results for linear programming method ###

C.labels <- 'would choose\n' %.% med_labels
names(C.labels) <- c('1', '2', '3')

treatment.labels <- c('13' = "a = pro-attitudinal\na' = entertainment",
                      '23' = "a = counter-attitudinal\na' = entertainment",
                      '12' = "a = pro-attitudinal\na' = counter-attitudinal")

acte.binary$sensitivity$effects$treatment <-
  factor(treatment.labels[acte.binary$sensitivity$effects$treat %.%
                            acte.binary$sensitivity$effects$ctrl],
         levels = treatment.labels)

## naive estimates (should be about the same as normal posterior)
acte.binary$naive.ci <- as.data.frame(acte.binary$naive$effects)
acte.binary$naive.ci$min <-
  acte.binary$posterior$naive_effects.ci$min_cilo
acte.binary$naive.ci$max <-
  acte.binary$posterior$naive_effects.ci$max_cihi
acte.binary$naive.ci$point <- NULL
acte.binary$naive.ci$treatment <-
  factor(treatment.labels[acte.binary$naive.ci$treat %.%
                            acte.binary$naive.ci$ctrl],
         levels = treatment.labels)

## bounds on treatment effects
acte.binary$effects.ci <- acte.binary$posterior$effects.ci
colnames(acte.binary$effects.ci) <-
  gsub('_ci(lo|hi)', '', colnames(acte.binary$effects.ci))
acte.binary$effects.ci$treatment <-
  factor(treatment.labels[acte.binary$effects.ci$treat %.%
                            acte.binary$effects.ci$ctrl],
         levels = treatment.labels)

## sensitivity analysis for treatment effects
acte.binary$sensitivity$effects.ci <-
  acte.binary$posterior$sens_effects.ci
colnames(acte.binary$sensitivity$effects.ci) <-
  gsub('_ci(lo|hi)', '', colnames(acte.binary$sensitivity$effects.ci))
acte.binary$sensitivity$effects.ci$treatment <-
  factor(treatment.labels[acte.binary$sensitivity$effects.ci$treat %.%
                            acte.binary$sensitivity$effects.ci$ctrl],
         levels = treatment.labels)

## naive estimates assume rho == 0
##   and bounds equivalent to rho == Inf
##   (for plotting purposes, place bounds at maximum rho)
acte.binary$naive$effects$rho <- 0
acte.binary$naive$effects$treatment <-
  factor(treatment.labels[acte.binary$naive$effects$treat %.%
                            acte.binary$naive$effects$ctrl],
         levels = treatment.labels)

acte.binary$naive.ci$rho <- 0
acte.binary$naive.ci$treatment <-
  factor(treatment.labels[acte.binary$naive.ci$treat %.%
                            acte.binary$naive.ci$ctrl],
         levels = treatment.labels)

acte.binary$effects$rho <-
  max(acte.binary$sensitivity$effects$rho[
    is.finite(acte.binary$sensitivity$effects$rho)
  ])
acte.binary$effects$treatment <-
  factor(treatment.labels[acte.binary$effects$treat %.%
                            acte.binary$effects$ctrl],
         levels = treatment.labels)

acte.binary$effects.ci$rho <-
  max(acte.binary$sensitivity$effects$rho[
    is.finite(acte.binary$sensitivity$effects$rho)
  ])
acte.binary$effects.ci$treatment <-
  factor(treatment.labels[acte.binary$effects.ci$treat %.%
                            acte.binary$effects.ci$ctrl])

pdf('./results/figA1_binary_sensitivity.pdf', 8, 8)
ggplot(acte.binary$sensitivity$effects[
  is.finite(acte.binary$sensitivity$effects$rho),], aes(x=rho)) +
  geom_ribbon(aes(ymin=min, ymax=max),
              data=acte.binary$sensitivity$effects.ci[
                is.finite(acte.binary$sensitivity$effects.ci$rho),],
              fill='#80808080') +
  geom_ribbon(aes(ymin=min, ymax=max),
              fill='#40404080') +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'naive'),
                data = acte.binary$naive.ci,
                width = .02, size = .5, linetype = 'dashed') +
  geom_point(aes(y = naive, colour = 'naive'),
             data = acte.binary$naive$effects) +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'bounds'),
                data = acte.binary$effects.ci,
                width = .015, size = .5) +
  geom_errorbar(aes(ymin = min, ymax = max, colour = 'bounds'),
                data = acte.binary$effects,
                width = .01, size = 1) +
  facet_grid(treatment ~ C, labeller = labeller(C = C.labels)) +
  scale_color_manual(values = c(naive = blue_medium, bounds = red_medium),
                     guide = FALSE) +
  xlab(expression('Sensitivity parameter ' * (rho))) +
  ylab(expression('Expected change in probability of discussing, ' * tau("a, a'| c"))) +
  scale_x_continuous(breaks = seq(0, .25, .1)) +
  scale_y_continuous(breaks=seq(-.5,.5,.25), limits=c(-.65,.65)) +
  ggtitle('Sensitivity Analysis for Discussing Story with Friends (binary)\n') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
