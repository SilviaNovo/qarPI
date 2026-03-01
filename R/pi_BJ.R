#' Prediction intervals based on the Box-Jenkins method
#'
#' Implements classical Box-Jenkins Gaussian prediction intervals for an AR(p) model
#' fitted by ordinary least squares.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param alpha Significance level. 1-\code{alpha} is the nominal coverage level. Default is 0.05.
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
#' This function implements the classical Box-Jenkins prediction intervals. The method fits an AR(p)
#' model to the observed time series via ordinary least squares (OLS) and derives analytical
#' prediction intervals under the assumption of Gaussian innovations.
#'
#' @references
#' Box, G. E. and Jenkins, G. M. (1976). Time Series Analysis: Forecasting and Control. Holden–Day,
#' San Francisco.
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
#' out <- pi_BJ(series,h=2)
#'
#' out$lpi
#' out$upi
#' out$len
#' @export
pi_BJ <- function(series, p = 1, h = 3, alpha = 0.05) {

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

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) stop("`alpha` must be a single number.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  if (n <= p) stop("`series` length must be greater than `p`.", call. = FALSE)

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))
  y <- series[(p + 1):n]

  fit <- stats::lm(y ~ xp)

  # Point forecasts
  x.fut <- numeric(p + h)
  x.fut[1:p] <- series[(n - p + 1):n]
  for (i in (p + 1):(p + h)) {
    x.fut[i] <- fit$coef[1] + fit$coef[-1] %*% x.fut[(i - 1):(i - p)]
  }

  # --------------------
  # MA(inf) representation
  # --------------------
  MAcoef <- function(phi, h) {
    psi <- numeric(h)
    psi[1] <- 1
    p <- length(phi)

    if (h > 1) {
      for (k in 2:h) {
        sum_val <- 0
        for (j in 1:min(p, k - 1)) {
          sum_val <- sum_val + phi[j] * psi[k - j]
        }
        psi[k] <- sum_val
      }
    }
    psi
  }

  coef <- fit$coef[-1]
  ma <- MAcoef(phi = coef, h = h)

  # --------------------
  # Box-Jenkins prediction interval
  # --------------------
  lpi.bj <- numeric(h)
  upi.bj <- numeric(h)

  sigma2 <- (n - 1) / n * stats::var(fit$residuals)
  z <- stats::qnorm(1 - alpha / 2)

  for (i in 1:h) {
    se_i <- sqrt(sigma2 * sum(ma[1:i]^2))
    lpi.bj[i] <- x.fut[i + p] - z * se_i
    upi.bj[i] <- x.fut[i + p] + z * se_i
  }

  l.bj <- upi.bj - lpi.bj

  list(pfor = x.fut[(p + 1):(p + h)], lpi = lpi.bj, upi = upi.bj, len = l.bj)
}

