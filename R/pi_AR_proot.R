#' Prediction intervals based on AR-proot algorithm
#'
#' Computes a bootstrap predictive-root-based prediction interval using the AR-proot algorithm for
#' an AR(p) model over a forecast horizon \code{h}.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param B Number of bootstrap replicates. Default is 1000.
#' @param alpha Significance level. 1-\code{alpha} is the nominal coverage level. Default is 0.05.
#' @param tau Quantile level used in estimation. Default is 0.5.
#'
#' @return A list with elements:
#' \describe{
#'   \item{pfor}{Numeric vector of point forecasts (length h).}
#'   \item{lpi}{Numeric vector of lower bounds (length \code{h}).}
#'   \item{upi}{Numeric vector of upper bounds (length \code{h}).}
#'   \item{len}{Numeric vector of interval lengths (length \code{h}).}
#' }
#'
#' @details
#' This function implements the AR-proot algorithm described in
#' Novo and Sánchez-Sellero (2025). Predictive residuals are first
#' obtained using a leave-one-out quantile autoregressive fit.
#' To account for the estimation uncertainty of the autoregressive
#' coefficients, bootstrap coefficient estimates are generated
#' through a multiplier (random weights) bootstrap scheme.
#' Conditional on these bootstrap coefficients, future bootstrap
#' predictions and bootstrap observations are constructed, and
#' predictive root replicates are defined as the difference between
#' bootstrap observations and their corresponding bootstrap predictions.
#' Equal-tailed prediction intervals are obtained from the empirical
#' quantiles of the bootstrap predictive roots.
#'
#' @examples
#' set.seed(123)
#'
#' # Simulation parameters
#' burn_in <- 300
#' n <- 25
#' m <- burn_in + n
#' coeff <- 0.6
#'
#' # Simulate AR(1) process
#' e <- rnorm(m)
#' x <- numeric(m)
#' x[1] <- rnorm(1)
#'
#' for (t in 2:m) {
#'   x[t] <- coeff * x[t - 1] + e[t]
#' }
#'
#' series <- ts(x[(burn_in + 1):m])
#'
#' # Compute prediction interval
#' out <- pi_AR_proot(series)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' @references
#' Novo, S., & Sánchez-Sellero, C. (2025).
#' \emph{Prediction intervals for quantile autoregression}.
#' arXiv:2512.22018. \url{https://arxiv.org/abs/2512.22018}
#'
#'
#' @export
pi_AR_proot<-function(series, p=1, h=3, B=1000,alpha=0.05,tau=0.5){
  # --------------------
  # Argument checks
  # --------------------
  if (missing(series)) stop("`series` is required.", call. = FALSE)
  if (!is.numeric(series) && !stats::is.ts(series)) {
    stop("`series` must be a numeric vector or a `ts` object.", call. = FALSE)
  }
  series <- as.numeric(series)
  if (anyNA(series)) stop("`series` contains NA values.", call. = FALSE)

  n <- length(series)
  if (n < 2L) stop("`series` must have length >= 2.", call. = FALSE)

  if (!is.numeric(p) || length(p) != 1L || is.na(p)) stop("`p` must be a single number.", call. = FALSE)
  if (p != as.integer(p) || p < 1L) stop("`p` must be a positive integer.", call. = FALSE)

  if (!is.numeric(h) || length(h) != 1L || is.na(h)) stop("`h` must be a single number.", call. = FALSE)
  if (h != as.integer(h) || h < 1L) stop("`h` must be a positive integer.", call. = FALSE)

  if (!is.numeric(B) || length(B) != 1L || is.na(B)) stop("`B` must be a single number.", call. = FALSE)
  if (B != as.integer(B) || B < 2L) stop("`B` must be an integer >= 2.", call. = FALSE)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) stop("`alpha` must be a single number.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau)) stop("`tau` must be a single number.", call. = FALSE)
  if (tau <= 0 || tau >= 1) stop("`tau` must be in (0, 1).", call. = FALSE)

  if (n <= p) stop("`series` length must be greater than `p`.", call. = FALSE)

  #Series tranformation in response+predictor format:
  #--------------------------------------------------

  n<-length(series)
  if(p==1) xp <- sapply(1:(n-p), function(x) series[(x+p-1):x])
  if(p>1) xp <- t(sapply(1:(n-p), function(x) series[(x+p-1):x]))
  y<-series[(p+1):n]

  x.p<-as.matrix(xp)

  #Predictive residuals
  #--------------------

  resid <- rep(0, (n - p))

  for (i in 1:(n - p)) {
    xtemp <- x.p[-i, , drop = FALSE]
    ytemp <- y[-i]
    tempdata <- data.frame(ytemp = ytemp, xtemp = xtemp)
    colnames(tempdata) <- c("ytemp", paste0("xtemp.", 1:p))
    tempfit <- quantreg::rq(ytemp ~ .,tau=tau, data = tempdata)
    newdata <- as.data.frame(as.list(x.p[i, ]))
    colnames(newdata) <- paste0("xtemp.", 1:p)
    resid[i] <- y[i] - stats::predict(tempfit, newdata = newdata)
  }
  res.p<-resid-stats::quantile(resid,tau)


  #Quantile model for tau:
  #----------------------------------------

  m.quantile<-quantreg::rq(y~xp,tau=tau)

  #Resampling:
  #-----------
  q.taup<-matrix(0,ncol=B,nrow=h)
  for(b in 1:B){
    set.seed(b)

    #Resampling residuals from their empiric distribution:
    #-----------------------------------------------------

    r.qar<-sample(res.p,h,replace=TRUE)

    #Boostrap future values. Boostrap version of the coefficients: Random weights bootstrap:
    #---------------------------------------------------------------------------------------
    w<-stats::rexp(n-p)
    m.quantile.star<-quantreg::rq(y~xp,tau=tau,weights=w)
    x.fut.star<-numeric(p+h)
    x.fut.star[1:p]<-series[(n-p+1):n]
    x.fut.star.hat<-numeric(p+h)
    x.fut.star.hat[1:p]<-series[(n-p+1):n]
    for (t in 1:h){
      x.fut.star.hat[t+p]<-m.quantile.star$coef[1]+m.quantile.star$coef[-1]%*%x.fut.star.hat[(p+t-1):t]
      x.fut.star[t+p]<-m.quantile$coef[1]+m.quantile$coef[-1]%*%x.fut.star[(p+t-1):t]+r.qar[t]
    }#end for t

    dif<-x.fut.star-x.fut.star.hat
    q.taup[,b]<-dif[(p+1):(h+p)]

  }#end for b

  #Prediction interval:
  #--------------------

  x.hat.fut<-numeric(p+h)
  x.hat.fut<-series[(n-p+1):n]
  for (t in (p+1):(p+h)) {x.hat.fut[t]=m.quantile$coef[1]+ m.quantile$coef[-1]%*%x.hat.fut[(t-1):(t-p)]}
  x.hat.f<-x.hat.fut[(p+1):(p+h)]

  q.alpha<-apply(q.taup,1,sort)[max(1L, floor(alpha/2*B)),]
  q.1alpha<-apply(q.taup,1,sort)[min(B, floor((1-alpha/2)*B)),]

  lpi.taup<-x.hat.f+q.alpha
  upi.taup<-x.hat.f+q.1alpha



  #Measures prediction interval:
  #-----------------------------
  l.taup=upi.taup-lpi.taup

  list(pfor=x.hat.f,lpi=lpi.taup, upi=upi.taup, len=l.taup)
}


