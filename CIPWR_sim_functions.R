inv.weibull = function(u, scale, shape){
  (-(log(u)/scale))^(1/shape)
}

simulate.data <- function(beta, alpha, gamma, n, d){
  
  nu <- 7
  lambda <- 0.01
  ntrt <- nrow(alpha)
  
  # simulate a data set with 3 treatment groups
  X1 <- rnorm(n,0,1)
  X2 <- rnorm(n,0,1)
  X3 <- rnorm(n,0,1)
  X4 <- rbinom(n,1,0.4)
  X5 <- runif(n,-2,2)
  
  X.all <- as.matrix(cbind(X1,X2,X3,X4,X5)) # all covariates
  X <- as.matrix(cbind(1, X1, X2, X3))
  V <- as.matrix(cbind(1, X1, X2, X4))
  W <- as.matrix(cbind(1, X1, X2, X5))
  
  ### generate treatment
  pi <- exp(V%*%t(alpha))/rowSums(exp(V%*%t(alpha)))
  Z <- Hmisc::rMultinom(pi, 1)
  Z <- as.factor(Z)
  Z.dummy <- t(sapply(Z, function(x) x==1:ntrt)*1)
  
  # Censoring time
  u <- runif(n,0,1)
  C.potential <- inv.weibull(u=u, scale=lambda*exp(W %*% t(gamma)), shape=nu)
  C <- rowSums(C.potential*Z.dummy)
  
  ### outcome
  eps <- matrix(rlogis(n*ntrt, 0, 7), nrow=n)
  T.potential <- X %*% t(beta) + eps
  T_ <- rowSums(T.potential*Z.dummy)
  Y.potential <- (T.potential>d)*1
  Y <- rowSums(Y.potential*Z.dummy)
  
  eta <- pmin(T_, C)
  DeltaT <- (T_<=C)*1
  DeltaC <- (C<=T_)*1
  R <- 1-(1-DeltaT)*(eta<d)
  Y[R==0] <- NA
  
  # for sufficiently large sample size
  mu.true <- 1-apply(Y.potential,2,mean)
  tau.true <- outer(mu.true, mu.true, "-")
  tau.true <- tau.true[lower.tri(tau.true)]
  
  return(list(X.all=X.all, Z=Z, Y=Y, eta=eta, DeltaT=DeltaT, DeltaC=DeltaC, C=C, C.potential=C.potential, T.potential=T.potential,
              R=R, mu.true=mu.true, ATE.true=tau.true))
}

simulate.data.randomC <- function(beta, alpha, gamma, n, d){
  
  nu <- 7
  lambda <- 0.01
  ntrt <- nrow(alpha)
  
  # simulate a data set with 3 treatment groups
  X1 <- rnorm(n,0,1)
  X2 <- rnorm(n,0,1)
  X3 <- rnorm(n,0,1)
  X4 <- rbinom(n,1,0.4)
  X5 <- runif(n,-2,2)
  
  X.all <- as.matrix(cbind(X1,X2,X3,X4,X5)) # all covariates
  X <- as.matrix(cbind(1, X1, X2, X3))
  V <- as.matrix(cbind(1, X1, X2, X4))
  W <- as.matrix(cbind(1, X1, X2, X5))
  
  ### generate treatment
  pi <- exp(V%*%t(alpha))/rowSums(exp(V%*%t(alpha)))
  Z <- Hmisc::rMultinom(pi, 1)
  Z <- as.factor(Z)
  Z.dummy <- t(sapply(Z, function(x) x==1:ntrt)*1)
  
  # Censoring time
  
  C1 <- runif(n, 20, 380)
  C2 <- runif(n, 30, 400)
  C3 <- runif(n, 50, 300)
  C <- rowSums(cbind(C1,C2,C3)*Z.dummy)
  
  ### outcome
  eps <- matrix(rlogis(n*ntrt, 0, 7), nrow=n)
  T.potential <- X %*% t(beta) + eps
  T_ <- rowSums(T.potential*Z.dummy)
  Y.potential <- (T.potential>d)*1
  Y <- rowSums(Y.potential*Z.dummy)
  
  eta <- pmin(T_, C)
  DeltaT <- (T_<=C)*1
  DeltaC <- (C<=T_)*1
  R <- 1-(1-DeltaT)*(eta<d)
  Y[R==0] <- NA
  
  # for sufficiently large sample size
  mu.true <- 1-apply(Y.potential,2,mean)
  tau.true <- outer(mu.true, mu.true, "-")
  tau.true <- tau.true[lower.tri(tau.true)]
  
  return(list(X.all=X.all, Z=Z, Y=Y, eta=eta, DeltaT=DeltaT, DeltaC=DeltaC, C=C, T.potential=T.potential,
              R=R, mu.true=mu.true, ATE.true=tau.true))
}

basehaz <- function(eta.sorted, eta.order, Delta.sorted.0, W, Z.sorted.j, gamma){
  ## Return the values of baseline hazard at each jump and their corresponding time and order
  DZ <- Delta.sorted.0 & Z.sorted.j
  t.jump <- eta.sorted[DZ] # ordered time where each jump occurs
  jump.idx <- eta.order[DZ] # indices of subjects who have jumps in the original vector of time (eta) 
  expterm <- exp(W %*% gamma)
  expterm <- expterm[eta.order]
  expterm <- expterm * as.numeric(Z.sorted.j)
  denom.tmp <- cumsum(expterm[length(expterm):1])[length(expterm):1] # denominator of the formula
  denom <- denom.tmp[DZ]
  return(list(t.jump=t.jump, basehaz=1/denom, idx=jump.idx))
}

basecumhaz <- function(t, basehaz.res){
  
  ### Estimate baseline cumulative hazard function (Breslow estimator) for group j
  # t: time at which baseline cumulative hazard is evaluated
  # basehaz.res: output from basehaz()
  
  sum(basehaz.res$basehaz[basehaz.res$t.jump<=t])
  
}

std.logis <- function(Y, Z, X, R, wts=NULL){
  
  tryCatch ( 
   {
      # This function is to calculate the standardized mean of fitted y from a logistic regression model
      # X: matrix of predictors for the outcomde model (intercept included)
      # Y: outcome variable
      # Z: treatment variable
      # R: indicator of noncensored observations I(C>min(T,d))
      # wts: total weights of each subject, i.e. 1/{pr(Z=z|x)*pr(R=1|Z=z,x)}
      
      ### Value
      # mhat is a data frame of the predicted probability for each treatment group
      # mu is a vector of standardized mean of fitted y for each treatment group
      # tau is a vector of estimated ATE (length of tau equals number of treatment pairs)
      
      ntrt <- length(table(Z))
      
      if(!all(X[,1]==1)){X <- as.matrix(cbind(1, X))}
      
      mhat <- matrix(NA, nrow=length(Z), ncol=ntrt) # predicted probability
      coef <- matrix(NA, nrow=ntrt, ncol=ncol(X))
      for(j in 1:ntrt){
        I <- R*(Z==j) # indicator of observing the potential outcome Y(k)
        yfit <- glm(Y~X-1, weights=I*wts, family="quasibinomial")
        mhat[,j] <- predict(yfit, newdata=data.frame(X), type="response")
        coef[j,] <- coef(yfit)
      }
      
      mu <- 1-unname(apply(mhat,2,mean))
      tau <- outer(mu, mu, "-")
      tau <- tau[lower.tri(tau)]
      allcombn <- combn(ntrt, 2)
      names(tau) <- apply(allcombn,2, function(x) paste0("EY(",x[2],")-EY(", x[1],")"))
      
      return(list(mhat=mhat, ATE=tau, mu=mu, coef=coef))
    }, warning=function(e) return(list(mhat=NA, ATE=NA, mu=NA, coef=NA)) )
}

s0 <- function(u,eta,W,Z,j,gamma){
  # u: time at which s0 is evaluated
  # eta: observed time
  # W: predictor matrix for censoring model, no intercept
  # Z: treatment variable
  # gamma.j: coefficients of censoring model for group j
  
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  n <- nrow(W)
  
  A <- (Z==j)*(eta>=u) # at risk indicator
  s0.out <- sum(exp(W%*%gamma.j)[A==1])/n
  return(s0.out)
  
}

s1 <- function(u,eta,W,Z,j,gamma){
  # u: time at which s1 is evaluated
  # eta: observed time
  # W: predictor matrix for censoring model
  # Z: treatment variable
  # j: treatment j
  # gamma.j: coefficients of censoring model for group j
  
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  n <- nrow(W)
  
  A <- (Z==j)*(eta>=u) # at risk indicator
  s1.out <- apply((as.numeric(exp(W%*%gamma.j))*W)[A==1,drop=F,],2,sum)/n
  return(s1.out)
  
}

s2 <- function(u,eta,W,Z,j,gamma){
  # u: time at which s2 is evaluated
  # eta: observed time, min(T,C)
  # W: predictor matrix for censoring model
  # Z: treatment variable
  # j: treatment j
  # gamma.j: coefficients of censoring model for group j
  
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  n <- nrow(W)
  
  A <- (Z==j)*(eta>=u) # at risk indicator
  
  s2.out <- (t(W*as.numeric(A*exp(W%*%gamma.j))) %*% W)/n
  
  return(s2.out)
  
}

estK.ij <- function(eta, xi, W, Z, j, gamma, basehaz.res, s0.res, s1.res){
  
  # Return an n by p matrix, where n is the total number of subjects and p is the number of predictors
  # in the censoring model.
  # eta: min(T,C)
  # xi: min(T,C,d)
  
  s1.res <- as.matrix(s1.res)
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  
  jump.idx <- basehaz.res$idx
  t.jump <- basehaz.res$t.jump
  basehaz <- basehaz.res$basehaz
  
  intensity <- matrix(0, nrow=n, ncol=ncol(W))
  
  for(i in 1:n){ 
    
    # calculate the intensity process part
    
    sum.tmp <- rep(0, ncol(W))
    jump.index <- 1
    
    if(Z[i]==j){
      while(t.jump[jump.index]<=xi[i] & jump.index<=length(t.jump)){
        integrand <- W[i,] - s1.res[jump.idx[jump.index],]/s0.res[jump.idx[jump.index]]
        sum.tmp <- sum.tmp + integrand*basehaz[jump.index]
        jump.index <- jump.index + 1
      }
    }
    intensity[i,] <- sum.tmp*as.numeric(exp(t(as.matrix(W[i,]))%*%as.matrix(gamma.j)))
  }
  
  return(intensity)
}

estP.j <- function(Z, j, X, K, R, Y, m, wts){
  n <- length(Y)
  D <- (Z==j)*1
  mj <- m[,j]
  v <- as.numeric(D*R*(Y-mj)*wts)
  v[is.na(v)] <- 0
  P.out <- (t(X*v) %*% K)/n
  
  return(P.out)
}

estQ.j <- function(Z, j, R, X, W, gamma, Y, m, wts){
  
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  
  n <- length(Y)
  D <- (Z==j)*1
  mj <- m[,j]
  v <- as.numeric(D*R*(Y-mj)*wts)
  v[is.na(v)] <- 0
  Q.out <- (t(X*v) %*% exp(W%*%gamma.j))/n
  return(Q.out)
}

estOmega.j <- function(eta, W, Z, j, gamma, basehaz.res, s0.res, s1.res){
  
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  
  # checked; results consistent with the output in SAS (model-based variance of coefficients)
  n <- nrow(W)
  t.jump <- basehaz.res$t.jump
  jump.idx <- basehaz.res$idx
  tmp.sum <- 0
  
  for(r in 1:length(t.jump)){
    s2hat <- s2(t.jump[r],eta,W,Z,j,gamma)
    s1hat <- as.matrix(s1.res[jump.idx[r],])
    s0hat <- as.numeric(s0.res[jump.idx[r]])
    tmp.sum <- tmp.sum + s2hat/s0hat - s1hat%*%t(s1hat)/(s0hat^2)
  }
  
  return(tmp.sum/n)
}

estRESOUT.ij <- function(X,Z,j,R,Y,m,wts){
  # 'residuals' of the outcome model
  mj <- m[,j]
  outcome.resid <- as.numeric((Z==j)*R*(Y-mj)*wts)*X
  outcome.resid[is.na(outcome.resid)] <- 0
  return(outcome.resid)
}

estRESPROPEN.il <- function(V,Z,l,gps){
  # gps is a matrix with (J-1) columns
  propen.resid <- as.numeric((Z==l)-gps[,l])*V
  return(propen.resid)
}

estU.ij <- function(eta, W, Z, j, gamma, basehaz.res, s0.res, s1.res){
  
  # estimate score residuals
  # checked; results consistent with SAS output
  
  s1.res <- as.matrix(s1.res)
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  
  n <- nrow(W)
  
  counting <- matrix(0, nrow=n, ncol=ncol(W))
  intensity <- matrix(0, nrow=n, ncol=ncol(W))
  
  jump.idx <- basehaz.res$idx
  basehaz <- basehaz.res$basehaz
  t.jump <- basehaz.res$t.jump
  counting[jump.idx,] <- W[jump.idx,] - s1.res[jump.idx,]/s0.res[jump.idx]
  
  for(i in 1:n){ 
    
    # calculate the intensity process part
    
    sum.tmp <- rep(0, ncol(W))
    jump.index <- 1
    
    if(Z[i]==j){
      while(t.jump[jump.index]<=eta[i] & jump.index<=length(t.jump)){
        integrand <- W[i,] - s1.res[jump.idx[jump.index],]/s0.res[jump.idx[jump.index]]
        sum.tmp <- sum.tmp + integrand*basehaz[jump.index]
        jump.index <- jump.index + 1
      }
    }
    intensity[i,] <- sum.tmp*as.numeric(exp(t(as.matrix(W[i,]))%*%as.matrix(gamma.j)))
  }
  
  score.resid <- counting - intensity
  return(score.resid)
  
}

estMART.ij <- function(eta, xi, W, Z, j, gamma, basehaz.res, s0.res){
  
  if(!is.matrix(W)){W <- as.matrix(W)}
  gamma.j <- gamma[j,]
  if(!is.matrix(gamma.j)){gamma.j <- as.matrix(gamma.j)}
  
  n <- nrow(W)
  
  # counting process
  counting <- rep(0, n)
  jump.idx <- basehaz.res$idx
  basehaz <- basehaz.res$basehaz
  t.jump <- basehaz.res$t.jump
  counting[jump.idx] <- 1/s0.res[jump.idx]
  counting[eta>xi] <- 0
  
  # intensity process
  integrand <- 1/s0.res[jump.idx]
  summation <- cumsum(integrand*basehaz)
  sum.idx <- sapply(xi, function(x) sum(t.jump<=x))
  sum.idx[sum.idx==0] <- NA
  summation.tmp <- summation[sum.idx]
  summation.tmp[is.na(summation.tmp)] <- 0
  summation.tmp[Z!=j] <- 0
  intensity <- as.numeric(summation.tmp*exp(W%*%gamma.j))
  
  mart <- counting - intensity
  return(mart)
}

estA.j <- function(X, j, m){
  # return a 1xp vector, where p is the number of columns in X
  # m: predicted probabilities obtained from outcome model
  mj <- m[,j]
  n <- nrow(X)
  A.out <- (t(mj*(1-mj)) %*% X)/n
  return(A.out)
}

estB.j <- function(X,Z,j,R,m,wts){
  mj <- m[,j]
  n <- nrow(X)
  v <- as.numeric((Z==j)*R*mj*(1-mj)*wts)
  B.out <- (t(X*v)%*%X)/n
  return(B.out)
}

estH.l <- function(V,l,gps){
  # gps: a matrix with J-1 columns, where J is the total number of treatment groups
  n <- nrow(V)
  v <- as.numeric(gps[,l]*(1-gps[,l]))
  H.out <- (t(V*v) %*% V)/n
  return(H.out)
}

estF.jl <- function(Z, j, l, Y, X, V, R, alpha, m, gps, surv){
  # gps: a matrix with (J-1) columns
  ntrt <- length(table(Z))
  mj <- m[,j]
  if(j==l){
    a <- 1-1/gps[,l]
  }else if(j==ntrt){
    denom <- 1
    a <- exp(V%*%alpha[l,])/denom
  } else if(j<ntrt){
    denom <- exp(V%*%alpha[j,])
    a <- exp(V%*%alpha[l,])/denom
  }
  surv.inv <- 1/surv
  surv.inv[which(!is.finite(surv.inv))] <- 0
  b <- as.numeric((Z==j)*R*(Y-mj)*a*surv.inv[,j])
  b[is.na(b)] <- 0
  F.out <- t(X*b) %*% V/n
  return(F.out)
}

estIPWR <- function(Y, Z, R, X, V, W, pi, surv, eta, Delta, gamma, d){
  
  # Y: outcome variable
  # Z: treatment variable
  # X: predictors for outcome model, intercept included
  # V: predictors for propensity score model, intercept included
  # W: predictors for outcome model, no intercept
  # pi: matrix of propensity scores, i.e. pr(Z=z|x)
  # surv: matrix of survival probability, i.e. pr(R=1|Z=z,x)
  # eta: min(T,C)
  # Delta: I(C<=T)
  # gamma: matrix of coefficients of the censoring model
  # d: prespecified cut-off time point
  
  tryCatch({
    Z <- factor(Z, levels=sort(unique(Z)))
    Z <- relevel(Z, ref=max(levels(Z)))
    
    if(!is.matrix(X)){X <- as.matrix(X)}
    if(!is.matrix(V)){V <- as.matrix(V)}
    if(!is.matrix(W)){W <- as.matrix(W)}
    
    n <- length(Y) # sample size
    
    ntrt <- length(table(Z))
    npair <- choose(ntrt, 2)
    
    xi <- pmin(eta,d)
    
    ### Obtain coefficients (alpha) of the propensity score model
    gps.fit <- nnet::multinom(Z~V-1, trace=F)
    gps <- fitted(gps.fit) # note the order of columns
    gps <- gps[,-1]
    alphahat <- coef(gps.fit)
    
    # Predicted probability P(Y=1|Z=z,x) for z=1,...,J
    Z.dummy <- t(sapply(Z, function(x) x==1:ntrt)*1)
    wts <- 1/rowSums(pi*surv*Z.dummy)
    stand.res <- std.logis(Y=Y, Z=Z, X=X, R=R, wts=wts)
    
    mhat <- stand.res$mhat
    tau <- stand.res$ATE
    mu <- stand.res$mu
    
    ###
    ### Treatment specific influence functions
    ###
    
    phi <- matrix(NA, nrow=n, ncol=ntrt)
    
    eta.sorted <- sort(eta)
    eta.order <- order(eta)
    Delta.sorted <- Delta[eta.order]
    Z.sorted <- Z[eta.order]
    
    for(j in 1:ntrt){
      s0.eta <- sapply(eta, function(x) s0(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma))
      s0.xi <- sapply(xi, function(x) s0(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma))
      if(dim(W)[2]>1){
        s1.eta <- t(sapply(eta, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
        s1.xi <- t(sapply(xi, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
      } else {
        s1.eta <- as.matrix(sapply(eta, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
        s1.xi <- as.matrix(sapply(xi, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
      }
      
      A <- estA.j(X=X, j=j, m=mhat)
      B <- estB.j(X=X, Z=Z, j=j, R=R, m=mhat, wts=wts)
      
      ## Residual terms
      resid1 <- mhat[,j]-mean(mhat[,j])
      
      RESOUT <- estRESOUT.ij(X, Z, j, R, Y, mhat, wts)
      resid2 <- RESOUT %*% MASS::ginv(B) %*% t(A)
      
      resid3 <- 0
      for(l in 1:(ntrt-1)){
        RESPROPEN <- estRESPROPEN.il(V=V,Z=Z,l=l,gps=gps)
        H <- estH.l(V=V,l=l,gps=gps)
        F.l <- estF.jl(Z=Z, j=j, l=l, Y=Y, X=X, V=V, R=R, alpha=alphahat, m=mhat, gps=gps, surv=surv)
        resid3 <- resid3 + RESPROPEN %*% MASS::ginv(H) %*% t(F.l) %*% MASS::ginv(B) %*% t(A)
      }
      
      basehaz.res <- basehaz(eta.sorted, eta.order, Delta.sorted==1, W, Z.sorted==j, gamma[j,])
      
      K <- estK.ij(eta, xi, W, Z, j, gamma, basehaz.res, s0.xi, s1.xi)
      P <- estP.j(Z=Z, j=j, X=X, K=K, R=R, Y=Y, m=mhat, wts=wts)
      Omega <- estOmega.j(eta, W, Z, j, gamma, basehaz.res, s0.eta, s1.eta)
      U <- estU.ij(eta, W, Z, j, gamma, basehaz.res, s0.eta, s1.eta)
      
      resid4 <- U %*% MASS::ginv(Omega) %*% t(P) %*% MASS::ginv(B) %*% t(A)
      
      Q <- estQ.j(Z, j, R, X, W, gamma, Y, mhat, wts)
      MART <- estMART.ij(eta, xi, W, Z, j, gamma, basehaz.res, s0.xi)
      resid5 <- as.numeric(t(Q) %*% MASS::ginv(B) %*% t(A))*MART
      
      ## Residual terms
      
      phi[,j] <- resid1 + resid2 + resid3 + resid4 + resid5
      
    }
    
    allcombn <- combn(ntrt, 2)
    
    var.tau <- NULL
    var.mu <- apply(phi,2,function(x) sum(x^2)/(n^2))
    
    for(kk in 1:npair){
      t1 <- allcombn[1,kk]; t2 <- allcombn[2,kk]
      var.tau[kk] <- sum((phi[,t1] - phi[,t2])^2)/(n^2)
    }
    
    mu.upper <- mu + qnorm(0.975)*sqrt(var.mu)
    mu.lower <- mu - qnorm(0.975)*sqrt(var.mu)
    
    tau.upper <- tau + qnorm(0.975)*sqrt(var.tau)
    tau.lower <- tau - qnorm(0.975)*sqrt(var.tau)
    
    
    return(list(mu=mu, ATE=tau, var.ATE=var.tau, var.mu=var.mu,
                ATE.upper=tau.upper, ATE.lower=tau.lower,
                mu.upper=mu.upper, mu.lower=mu.lower))
  }, error=function(e) return(list(mu=NA, ATE=NA, var.ATE=NA, var.mu=NA,
                                   ATE.upper=NA, ATE.lower=NA,
                                   mu.upper=NA, mu.lower=NA)) )
}

estIPWR.C <- function(Y, Z, R, X, V, W, pi, surv, eta, Delta, gamma, d){
  
  # Y: outcome variable
  # Z: treatment variable
  # X: predictors for outcome model, intercept included
  # V: predictors for propensity score model, intercept included
  # W: predictors for outcome model, no intercept
  # pi: matrix of propensity scores, i.e. pr(Z=z|x)
  # surv: matrix of survival probability, i.e. pr(R=1|Z=z,x)
  # eta: min(T,C)
  # Delta: I(C<=T)
  # gamma: matrix of coefficients of the censoring model
  # d: prespecified cut-off time point
  
  tryCatch({
  Z <- factor(Z, levels=sort(unique(Z)))
  Z <- relevel(Z, ref=max(levels(Z)))
  
  if(!is.matrix(X)){X <- as.matrix(X)}
  if(!is.matrix(V)){V <- as.matrix(V)}
  if(!is.matrix(W)){W <- as.matrix(W)}
  
  n <- length(Y) # sample size
  
  ntrt <- length(table(Z))
  npair <- choose(ntrt, 2)
  
  xi <- pmin(eta,d)
  
  ### Obtain coefficients (alpha) of the propensity score model
  gps.fit <- nnet::multinom(Z~V-1, trace=F)
  gps <- fitted(gps.fit) # note the order of columns
  gps <- gps[,-1]
  alphahat <- coef(gps.fit)
  
  # Predicted probability P(Y=1|Z=z,x) for z=1,...,J
  Z.dummy <- t(sapply(Z, function(x) x==1:ntrt)*1)
  wts <- 1/rowSums(pi*surv*Z.dummy)
  stand.res <- std.logis(Y=Y, Z=Z, X=X, R=R, wts=wts)
  
  mhat <- stand.res$mhat
  tau <- stand.res$ATE
  mu <- stand.res$mu
  
  ###
  ### Treatment specific influence functions
  ###
  
  phi <- matrix(NA, nrow=n, ncol=ntrt)
  
  eta.sorted <- sort(eta)
  eta.order <- order(eta)
  Delta.sorted <- Delta[eta.order]
  Z.sorted <- Z[eta.order]
  
  for(j in 1:ntrt){
    s0.eta <- sapply(eta, function(x) s0(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma))
    s0.xi <- sapply(xi, function(x) s0(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma))
    if(dim(W)[2]>1){
      s1.eta <- t(sapply(eta, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
      s1.xi <- t(sapply(xi, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
    } else {
      s1.eta <- as.matrix(sapply(eta, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
      s1.xi <- as.matrix(sapply(xi, function(x) s1(x, eta=eta, W=W, Z=Z, j=j, gamma=gamma)))
    }
    
    A <- estA.j(X=X, j=j, m=mhat)
    B <- estB.j(X=X, Z=Z, j=j, R=R, m=mhat, wts=wts)
    
    ## Residual terms
    resid1 <- mhat[,j]-mean(mhat[,j])
    
    RESOUT <- estRESOUT.ij(X, Z, j, R, Y, mhat, wts)
    resid2 <- RESOUT %*% MASS::ginv(B) %*% t(A)
    
    resid3 <- 0
    for(l in 1:(ntrt-1)){
      RESPROPEN <- estRESPROPEN.il(V=V,Z=Z,l=l,gps=gps)
      H <- estH.l(V=V,l=l,gps=gps)
      F.l <- estF.jl(Z=Z, j=j, l=l, Y=Y, X=X, V=V, R=R, alpha=alphahat, m=mhat, gps=gps, surv=surv)
      resid3 <- resid3 + RESPROPEN %*% MASS::ginv(H) %*% t(F.l) %*% MASS::ginv(B) %*% t(A)
    }
    
    basehaz.res <- basehaz(eta.sorted, eta.order, Delta.sorted==1, W, Z.sorted==j, gamma[j,])
    
    K <- estK.ij(eta, xi, W, Z, j, gamma, basehaz.res, s0.xi, s1.xi)
    P <- estP.j(Z=Z, j=j, X=X, K=K, R=R, Y=Y, m=mhat, wts=wts)
    Omega <- estOmega.j(eta, W, Z, j, gamma, basehaz.res, s0.eta, s1.eta)
    U <- estU.ij(eta, W, Z, j, gamma, basehaz.res, s0.eta, s1.eta)
    
    resid4 <- U %*% MASS::ginv(Omega) %*% t(P) %*% MASS::ginv(B) %*% t(A)
    
    Q <- estQ.j(Z, j, R, X, W, gamma, Y, mhat, wts)
    MART <- estMART.ij(eta, xi, W, Z, j, gamma, basehaz.res, s0.xi)
    resid5 <- as.numeric(t(Q) %*% MASS::ginv(B) %*% t(A))*MART
    
    ## Residual terms
    
    phi[,j] <- resid1 + resid2 + resid3 + resid4 + resid5
    
  }
  
  allcombn <- combn(ntrt, 2)
  
  var.tau <- NULL
  var.mu <- apply(phi,2,function(x) sum(x^2)/(n^2))
  
  for(kk in 1:npair){
    t1 <- allcombn[1,kk]; t2 <- allcombn[2,kk]
    var.tau[kk] <- sum((phi[,t1] - phi[,t2])^2)/(n^2)
  }
  
  mu.upper <- mu + qnorm(0.975)*sqrt(var.mu)
  mu.lower <- mu - qnorm(0.975)*sqrt(var.mu)
  
  tau.upper <- tau + qnorm(0.975)*sqrt(var.tau)
  tau.lower <- tau - qnorm(0.975)*sqrt(var.tau)
  
  
  return(list(mu=mu, ATE=tau, var.ATE=var.tau, var.mu=var.mu,
              ATE.upper=tau.upper, ATE.lower=tau.lower,
              mu.upper=mu.upper, mu.lower=mu.lower))
  }, error=function(e) return(list(mu=NA, ATE=NA, var.ATE=NA, var.mu=NA,
                                   ATE.upper=NA, ATE.lower=NA,
                                   mu.upper=NA, mu.lower=NA)) )
}


IPW <- function(Y, R, Z, wts){
  
  ntrt <- length(table(Z))
  mu <- sapply(1:ntrt, function(j) 1-weighted.mean(Y, R*(Z==j)*wts, na.rm=T))
  tau <- outer(mu, mu, "-")
  tau <- tau[lower.tri(tau)]
  return(list(mu=mu, ATE=tau))
}

KM <- function(eta, Delta, Z=NULL, j=NULL){
  # return KM estimates and the corresponding time
  if(!is.null(Z)){eta=eta[Z==j]; Delta=Delta[Z==j]}
  eta.sorted <- sort(eta)
  eta.order <- order(eta)
  Delta.sorted <- Delta[eta.order]
  n <- length(eta.sorted)
  atrisk <- n:1
  S <- cumprod((1-Delta.sorted/atrisk))
  #S.ori <- S
  #S.ori[eta.order] <- S
  return(list(surv=S, time=eta.sorted))
}

AIPW.Wang <- function(Z, Y, X, gps, R, surv){
  tryCatch ( 
    {
  ntrt <- length(table(Z))
  
  if(!all(X[,1]==1)){X <- as.matrix(cbind(1, X))}
  w <- R/surv
  mu <- NULL
  for(j in 1:ntrt){
    yfit <- glm(Y~X-1, weights=(Z==j)*w, family="quasibinomial")
    mhat <- predict(yfit, newdata=data.frame(X), type="response")
    x <- (Z==j)*Y/gps[,j] - ((Z==j)-gps[,j])*mhat/gps[,j]
    mu[j] <- weighted.mean(x,w)
  }
  mu <- 1-mu
  tau <- outer(mu, mu, "-")
  tau <- tau[lower.tri(tau)]
  allcombn <- combn(ntrt, 2)
  names(tau) <- apply(allcombn,2, function(x) paste0("EY(",x[2],")-EY(", x[1],")"))
  return(list(mu=mu, ATE=tau))
    }, error=function(e) return(list(mu=NA, ATE=NA)))
}


estG <- function(t, basehazT.res, basehazC.res, exptermT, exptermC, C, ps, trtind){
  
  n <- length(exptermT)
  idxC <- basehazC.res$idx # indices of subjects that are censored
  basehazC <- basehazC.res$basehaz
  tJumpC <- basehazC.res$t.jump
  
  basecumC <- sapply(tJumpC[tJumpC<=t], function(x) basecumhaz(x, basehazC.res))
  basecumT <- sapply(tJumpC[tJumpC<=t], function(x) basecumhaz(x, basehazT.res))
  
  # counting process
  G1 <- matrix(0, nrow=sum(tJumpC<=t), ncol=n)
  G1.rows <- 1:sum(tJumpC<=t); G1.cols <- idxC[tJumpC<=t]
  cumhazC <- basecumC * exptermC[idxC[tJumpC<=t]]
  cumhazT <- basecumT * exptermT[idxC[tJumpC<=t]]
  G1[cbind(G1.rows, G1.cols)] <- exp(cumhazC + cumhazT)
  G1 <- apply(G1, 2, cumsum)
  
  # intensity process
  
  integrand <- matrix(exp(exptermC %o% basecumC + exptermT %o% basecumT), nrow=n)
  A <- outer(C, tJumpC[tJumpC<=t], ">=")*trtind # at risk indicator
  jumpSum <- t(sweep(integrand*A, MARGIN=2, STATS=basehazC[tJumpC<=t], FUN='*'))
  intensity <- sweep(jumpSum, MARGIN=2, STATS=exptermC, FUN='*')
  G2 <- apply(intensity,2,cumsum)
  
  G <- sweep(G1-G2, MARGIN=2, STATS=trtind*1/ps, FUN='*')
  
  return(G)
}

surv.ZS <- function(t, basehazT.res, basehazC.res, exptermT, exptermC, eta, C, ps, trtind){
  tryCatch({
    n <- length(exptermT)
    idxC <- basehazC.res$idx
    basehazC <- basehazC.res$basehaz
    tJumpC <- basehazC.res$t.jump
    idxT <- basehazT.res$idx
    basehazT <- basehazT.res$basehaz
    tJumpT <- basehazT.res$t.jump
    
    # denominator
    w <- (trtind/ps)*1
    basecumC <- sapply(tJumpT[tJumpT<=t], function(x) basecumhaz(x, basehazC.res)) # baseline cum hazard of C evaluated at tJump.T
    basecumT <- sapply(tJumpT[tJumpT<=t], function(x) basecumhaz(x, basehazT.res))
    
    summand1 <- w*matrix(exp(exptermC %o% basecumC), nrow=n)
    A <- outer(eta, tJumpT[tJumpT<=t], ">=")*trtind # at risk indicator
    denom1 <- apply(summand1*A, 2, sum)
    
    G <- estG(t, basehazT.res, basehazC.res, exptermT, exptermC, C, ps, trtind)
    G <- rbind(0, G)
    summand2 <- matrix(exp(-(exptermT %o% basecumT)), nrow=n)
    denom2 <- apply(summand2, 2, sum)
    
    summand3 <- -w*matrix(exp(-(exptermT %o% basecumT)), nrow=n)
    denom3 <- apply(summand3, 2, sum)
    
    summand4 <- matrix(exp(-(exptermT %o% basecumT)), nrow=n)
    idx <- sapply(tJumpT[tJumpT<=t], function(x) sum(x>=tJumpC[tJumpC<=t]))
    GT <- t(G[idx+1,])
    GT[is.na(GT)] <- 0
    denom4 <- apply(summand4*GT, 2, sum)
    
    denom <- denom1 + denom2 + denom3 + denom4
    
    summand5 <- matrix(0, nrow=n, ncol=sum(tJumpT<=t))
    idxCols <- 1:sum(tJumpT<=t); idxRows <- idxT[tJumpT<=t]
    summand5[cbind(idxRows, idxCols)] <- w[idxT[tJumpT<=t]]*exp(basecumC * exptermC[idxT[tJumpT<=t]])
    nume1 <- apply(summand5, 2, sum)
    
    summand6 <- matrix(exp(-exptermT%o%basecumT), nrow=n) * as.vector(exptermT)
    summand6 <- sweep(summand6, MARGIN=2, STATS=basehazT[tJumpT<=t], FUN="*")
    nume2 <- apply(summand6, 2, sum)
    
    summand7 <- -matrix(exp(-exptermT%o%basecumT), nrow=n) * as.vector(exptermT*w)
    summand7 <- sweep(summand7, MARGIN=2, STATS=basehazT[tJumpT<=t], FUN="*")
    nume3 <- apply(summand7, 2, sum)
    
    summand8 <- GT * matrix(exp(-exptermT%o%basecumT), nrow=n) * as.vector(exptermT)
    summand8 <- sweep(summand8, MARGIN=2, STATS=basehazT[tJumpT<=t], FUN="*")
    nume4 <- apply(summand8, 2, sum)
    
    nume <- nume1 + nume2 + nume3 + nume4
    
    cumhaz.AIPW <- ifelse(sum(nume/denom)>=0, sum(nume/denom), 0)
    
    cumhaz.IPW <- ifelse(sum(nume1/denom1)>=0, sum(nume1/denom1), 0)
    
    return(list(cumhaz.IPW=cumhaz.IPW, surv.IPW=exp(-cumhaz.IPW), 
                cumhaz.AIPW=cumhaz.AIPW, surv.AIPW=exp(-cumhaz.AIPW)))
  }, error=function(err) return(list(cumhaz.IPW=NA, surv.IPW=NA, 
                                     cumhaz.AIPW=NA, surv.AIPW=NA)))
  
}

