# Workhorse function
lll <- function(Y, D, A, X1, X2, y.model, tol){

  ## Implement Long, Little and Lin's (2008) EM algorithm
  tlevels <- levels(factor(A))
  J <- length(tlevels)

  # Get intitial parameter values from naive fits
  data.c <- data.frame(C=factor(A), X2=X2, D=factor(D), w=1)
  data.c.mlogit <- mlogit.data(data.c, shape="wide", choice="C")
  fit.c <- mlogit(C ~ 0 | X2, data = data.c.mlogit[data.c.mlogit$D==0,], weights=w)
  C.init <- factor(sample(tlevels, length(Y), replace=TRUE))
  data.y <- data.frame(Y=Y, A=factor(A), X2=X2, D=factor(D), C=C.init)
  if(!is.null(X1)){
    data.y$X1 <- X1
  }
  data.y$C[data.y$D==0] <- data.y$A[data.y$D==0]
  if(!is.null(X1)){
    form.y <- Y ~ X1 + A*C
  } else {
    form.y <- Y ~ A*C
  }
  if(y.model == "normal"){
    fit.y <- lm(form.y, data = data.y)
  } else if (y.model == "logit"){
    fit.y <- glm(form.y, data = data.y, family=binomial(link="logit"))
  } else stop("incorrect y model")

  # EM loop
  maxdiff <- 100000000
  count <- 0
  while(maxdiff > tol & count < 100){

    # E step (calculate weights)
    pc.s <- predict(fit.c, newdata = data.c.mlogit[data.c.mlogit$D==1,])
    sum.fit.y <- summary(fit.y)
    data.y.aug <- fy.asc <- list(rep(NA, J))
    wnum <- matrix(NA, nrow=nrow(pc.s), ncol=J)
    for(j in 1:J){
      data.y.aug[[j]] <- subset(data.y, D==1)
      data.y.aug[[j]]$C <- tlevels[j]
      means.y.asc <- predict(fit.y, newdata = data.y.aug[[j]], type = "response")
      if(y.model == "normal"){
        fy.asc[[j]] <- dnorm(data.y.aug[[j]]$Y, mean=means.y.asc, sd=sum.fit.y$sigma)
      } else if (y.model == "logit"){
        fy.asc[[j]] <- dbinom(data.y.aug[[j]]$Y, 1, prob = means.y.asc)
      } else stop("incorrect y model")
      wnum[,j] <- fy.asc[[j]]*pc.s[,j]
    }
    w <- t(apply(wnum, 1, function(x) x/sum(x)))
    for(j in 1:J){
      data.y.aug[[j]] <- cbind(data.y.aug[[j]], w=w[,j])
    }

    # M step (fit weighted models)
    data.y.aug.stack <- do.call(rbind, data.y.aug)
    data.y.0 <- data.y[data.y$D==0,]
    data.y.0$w <- 1
    data.y.aug.stack <- rbind(data.y.aug.stack, data.y.0)

    data.c.aug.stack <- subset(data.y.aug.stack, select = c("C","X2","D","w"))
    data.c.aug.mlogit <- mlogit.data(data.c.aug.stack, shape = "wide", choice = "C")
    fit.c.new <- mlogit(C ~ 0 | X2, data = data.c.aug.mlogit, weights = w)
    if(y.model == "normal"){
      fit.y.new <- lm(form.y, data = data.y.aug.stack, weights = w)
    } else if (y.model == "logit"){
      fit.y.new <- glm(form.y, data = data.y.aug.stack, family=binomial(link="logit"), weights = w)
    } else stop("incorrect y model")

    # update parameters and assess convergence
    diff <- c(coef(fit.y.new) - coef(fit.y), coef(fit.c.new) - coef(fit.c))
    if(y.model == "normal"){
      diff <- c(diff, sigma(fit.y.new) - sigma(fit.y))
    }
    maxdiff <- max(abs(diff))
    fit.c <- fit.c.new
    fit.y <- fit.y.new
    cat("updated max difference = ", maxdiff, "\n")
    count <- count + 1
    if(count == 100) cat("convergence not achieved! \n")
  }

  pred <- as.data.frame(
    prop.table(table(X1 = data.y$X1, X2 = data.y$X2, C = data.y$C))
  )
  for (a in 1:J){
    pred$A = factor(a, 1:J)
    pred[,'A.' %.% a] <- predict(fit.y, pred)
  }
  pred$A <- NULL
  strata <- ddply(pred, 'C', function(pred.c){
    apply(pred.c[,'A.' %.% 1:J], 2, function(yhat.c){
      weighted.mean(yhat.c, pred.c$Freq)
    })
  })
  rownames(strata) <- 'C.' %.% strata$C
  strata$C <- NULL
  strata <- as.matrix(strata)
  return(list(fit.c = fit.c, fit.y = fit.y, strata = strata))

}

## # wrapper for boot
## lll.boot <- function(dd, ind, y.model, tol){
##   Y <- dd$Y[ind]
##   D <- dd$D[ind]
##   A <- dd$A[ind]
##   X1 <- dd$X1[ind]
##   X2 <- dd$X2[ind]
##   lllout <- lll(Y=Y, D=D, A=A, X1=X1, X2=X2, y.model=y.model, tol=tol)
##   alpha <- coef(lllout$fit.c)
##   alpha.se <- summary(lllout$fit.c)$CoefTable[,2]
##   beta <- coef(lllout$fit.y)
##   beta.se <- sqrt(diag(summary(lllout$fit.y)$cov.unscaled))
##   estimates <- c(beta, alpha, beta.se, alpha.se)
##   return(estimates)
## }

## ## Replicate LLL simulation
## lll.repl <- function(x, J, N, R, alpha, beta, tol, ncpus.boot){
##   X2 <- runif(N)
##   C.lin <- crossprod(rbind(rep(1,N), X2), alpha)
##   C <- rbinom(N, 1, plogis(C.lin))
##   X1 <- runif(N)
##   D <- rbinom(N, 1, .5)
##   A <- (D == 1)*rbinom(N, 1, .5) + (D == 0)*C
##   Y.lin <- crossprod(rbind(rep(1,N), X1, A, C, A*C), beta)
##   Y <- rbinom(N, 1, plogis(Y.lin))
##   dd <- data.frame(Y, D=D, A=A, X1=X1, X2=X2)

##   out.boot <- boot(dd, lll.boot, R=R, y.model="logit", tol=tol, parallel="multicore", ncpus=ncpus.boot)
##   est.coef <- out.boot$t0[1:(length(alpha)+length(beta))]
##   est.se.info <- out.boot$t0[(length(alpha)+length(beta)+1):(2*(length(alpha)+length(beta)))]
##   est.se.boot <- apply(out.boot$t[,1:(length(alpha)+length(beta))], 2, sd)
##   out.repl <- c(est.coef, est.se.info, est.se.boot)
##   return(out.repl)
## }

## J <- 2
## N <- 1000
## R <- 200
## alpha <- c(-1,2)
## beta <- c(-2, 2, 2, 2, -2)
## ncpus.boot <- 1
## nrepl <- 100
## ncpus.sf <- 20
## tol <- 0.001

## sfInit(parallel=TRUE, cpus=ncpus.sf)
## sfLibrary(boot)
## sfLibrary(mlogit)
## sfExport("J","N","R","alpha","beta","tol","ncpus.boot","lll","lll.boot","lll.repl")
## lll.repl.out <- sfSapply(1:nrepl, lll.repl, J=J, N=N, R=R, alpha=alpha, beta=beta, tol=tol, ncpus.boot=ncpus.boot)
## sfStop()

## lll.repl.tab <- matrix(NA, nrow=7, ncol=5,
##                        dimnames = list(c("Intercept","X1","A","C","A*C","Intercept","X2"),
##                                        c("True Value","Estimate","SE Info", "SE Boot", "Empirical SE")))
## lll.repl.tab[,1] <- c(beta,alpha)
## lll.repl.tab[,2] <- apply(lll.repl.out, 1, mean)[1:7]
## lll.repl.tab[,3] <- apply(lll.repl.out, 1, mean)[8:14]
## lll.repl.tab[,4] <- apply(lll.repl.out, 1, mean)[15:21]
## lll.repl.tab[,5] <- apply(lll.repl.out, 1, sd)[1:7]
## print(round(lll.repl.tab, 2))

## ## testing a setup similar to ours

## lll.ours <- function(x, J, N, R, alpha, beta, tol, ncpus.boot){
##   D <- rbinom(N, 1, .5)
##   S.dum <- rmultinom(N, 1, rep(1,J))
##   S <- factor(crossprod(S.dum, 1:J))
##   A.dum <- rmultinom(N, 1, rep(1,J))
##   A <- factor(crossprod(A.dum, 1:J))
##   C.lin <- crossprod(rbind(rep(1,N), S.dum[2:J,]), alpha)
##   C.mean <- exp(C.lin)/apply(exp(C.lin),1,sum)
##   C.dum <- apply(C.mean, 1, function(x) rmultinom(1, 1, x))
##   C <- factor(crossprod(C.dum, 1:J))
##   A[D==0] <- C[D==0]

##   X <- model.matrix(~ S + A*C)
##   Y.lin <- X %*% beta
##   Y.bin <- rbinom(N, 1, plogis(Y.lin))
##   Y.cont <- rnorm(N, Y.lin, 1)

##   dd <- data.frame(Y=Y.bin, D=factor(D), A=factor(A), X1=factor(S), X2=factor(S))

##   out.boot <- boot(dd, lll.boot, R=R, y.model="logit", tol=tol, parallel="multicore", ncpus=ncpus.boot)
##   est.coef <- out.boot$t0[1:(length(alpha)+length(beta)-J)]
##   est.se.info <- out.boot$t0[(length(alpha)+length(beta)+1-J):(2*(length(alpha)+length(beta)-J))]
##   est.se.boot <- apply(out.boot$t[,1:(length(alpha)+length(beta)-J)], 2, sd)
##   out.repl <- c(est.coef, est.se.info, est.se.boot)
##   return(out.repl)
## }

## J <- 3
## N <- 1000
## R <- 200
## alpha <- matrix(c(0,0,0, -0.5,2,1, -1,0.5,3), nrow=3, ncol=3)
## #beta <- c(1, 0, 0, 1.5, -1.5, -1, 1, 0.5, -0.5, -0.8, 0.8)  # true model for X1=NULL
## beta <- c(1, 2, -2, 1.5, -1.5, -1, 1, 0.5, -0.5, -0.8, 0.8)  # true model for X1!=NULL
## tol <- 0.01
## ncpus.boot <- 1
## nrepl <- 100
## ncpus.sf <- 20

## sfInit(parallel=TRUE, cpus=ncpus.sf)
## sfLibrary(boot)
## sfLibrary(mlogit)
## sfExport("J","N","R","alpha","beta","tol","ncpus.boot","lll","lll.boot","lll.ours")
## lll.ours.out <- sfSapply(1:nrepl, lll.ours, J=J, N=N, R=R, alpha=alpha, beta=beta, tol=tol, ncpus.boot=ncpus.boot)
## sfStop()

## lll.ours.tab <- matrix(NA, nrow=17, ncol=5,
##                        dimnames = list(c("Intercept","X12","X13","A2","A3","C2","C3","A2:C2","A3:C2","A2:C3","A3:C3",
##                                          "2:Intercept","3:Intercept","2:X22","3:X22","2:X23","3:X23"),
##                                        c("True Value","Estimate","SE Info", "SE Boot", "Empirical SE")))
## lll.ours.tab[,1] <- c(beta,t(alpha[,2:3]))
## lll.ours.tab[,2] <- apply(lll.ours.out, 1, mean)[1:17]
## lll.ours.tab[,3] <- apply(lll.ours.out, 1, mean)[18:34]
## lll.ours.tab[,4] <- apply(lll.ours.out, 1, mean)[35:51]
## lll.ours.tab[,5] <- apply(lll.ours.out, 1, sd)[1:17]
## print(round(lll.ours.tab, 2))

## # summary(mlogit(C ~ 0 | S, data = dd, shape = "wide"))
## # summary(lm(Y.cont ~ S + A*C, data = dd))
## # summary(glm(Y.bin ~ S + A*C, data = dd, family=binomial(link="logit")))
## #
## # out.bin <- lll(Y=Y.bin, D=factor(D), A=factor(A), X1=factor(S), X2=factor(S), y.model="logit", tol=0.001)
## # out.cont <- lll(Y=Y.cont, D=factor(D), A=factor(A), X1=factor(S), X2=factor(S), y.model="normal", tol=0.001)

