#' Prediction intervals based on the PRR algorithm
#'
#' Computes bootstrap percentile-based prediction intervals using the forward algorithm of
#' Pascual, Romo and Ruiz (2004) for an AR(p) model over a
#' forecast horizon \code{h}.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param B Number of bootstrap replicates. Default is 1000.
#' @param alpha Significance level. 1-\code{alpha} is the nominal coverage level. Default is 0.05.
#' @param method Estimation method. One among "OLS" and "LAD". Default is "LAD".
#'
#' @return A list with elements:
#' \describe{
#'   \item{bfor}{Numeric matrix of bootstrap forecasts with dimension h x B.}
#'   \item{lpi}{Numeric vector of lower bounds (length \code{h}).}
#'   \item{upi}{Numeric vector of upper bounds (length \code{h}).}
#'   \item{len}{Numeric vector of interval lengths (length \code{h}).}
#' }
#'
#' @details
#' This function implements the forward bootstrap algorithm described in
#' Pascual, Romo and Ruiz (2004).
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
#' out <- pi_PRR(series,h=4)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' @references
#' Pascual, L., Romo, J., and Ruiz, E. (2004). Bootstrap predictive inference for ARIMA processes.
#' \emph{Journal of Time Series Analysis}, 25(4):449–465.
#'
#' @export

pi_PRR <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05,
                   method = c("LAD", "OLS")) {

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
  if (B != as.integer(B) || B < 2L) stop("`B` must be a positive integer >= 2.", call. = FALSE)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) stop("`alpha` must be a single number.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  method <- match.arg(method)

  if (n <= p) stop("`series` length must be greater than `p`.", call. = FALSE)
  if (n - 2L * p <= 0L) stop("Need `n - 2*p > 0` for the residual rescaling used in this function.", call. = FALSE)

  if (method == "LAD" && !requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for method='LAD'. Install it with install.packages('quantreg').",
         call. = FALSE)
  }

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))
  y <- series[(p + 1):n]

  # --------------------
  # Estimation + residual rescaling
  # --------------------
  if (method == "LAD") {
    m.lad <- quantreg::rq(y ~ xp, tau = 0.5)
    r.lad <- sqrt((n - p) / (n - 2 * p)) * (m.lad$residuals - stats::median(m.lad$residuals))
  }

  if (method == "OLS") {
    m.ols <- stats::lm(y ~ xp)
    r.ols.res <- sqrt((n - p) / (n - 2 * p)) * (m.ols$residuals - base::mean(m.ols$residuals))
  }

  # --------------------
  # Resampling
  # --------------------
  p.prr <- matrix(0, ncol = B, nrow = h)

  for (b in 1:B) {
    set.seed(b)

    if (method == "OLS") {
      x.fut <- numeric(n + h)
      x.fut[1:p] <- series[1:p]

      res.ols <- sample(r.ols.res, n + h, replace = TRUE)

      # Forward bootstrap resampling (build bootstrap series up to n)
      for (i in (p + 1):n) {
        x.fut[i] <- m.ols$coef[1] + m.ols$coef[-1] %*% x.fut[(i - 1):(i - p)] + res.ols[i]
      }

      # Bootstrap version of the coefficients
      if (p == 1) xp.star <- sapply(1:(n - p), function(x) x.fut[(x + p - 1):x])
      if (p > 1)  xp.star <- t(sapply(1:(n - p), function(x) x.fut[(x + p - 1):x]))
      y.star <- x.fut[(p + 1):n]
      m.ols.star <- stats::lm(y.star ~ xp.star)

      # Bootstrap future values (use original last p values)
      x.fut[n:(n - p + 1)] <- series[n:(n - p + 1)]
      for (t in (n + 1):(n + h)) {
        x.fut[t] <- m.ols.star$coef[1] + m.ols.star$coef[-1] %*% x.fut[(t - 1):(t - p)] + res.ols[t]
      }
    }

    if (method == "LAD") {
      x.fut <- numeric(n + h)
      x.fut[1:p] <- series[1:p]

      res.lad <- sample(r.lad, n + h, replace = TRUE)

      # Forward bootstrap resampling (build bootstrap series up to n)
      for (i in (p + 1):n) {
        x.fut[i] <- m.lad$coef[1] + m.lad$coef[-1] %*% x.fut[(i - 1):(i - p)] + res.lad[i]
      }

      # Bootstrap version of the coefficients
      if (p == 1) xp.star <- sapply(1:(n - p), function(x) x.fut[(x + p - 1):x])
      if (p > 1)  xp.star <- t(sapply(1:(n - p), function(x) x.fut[(x + p - 1):x]))
      y.star <- x.fut[(p + 1):n]
      m.lad.star <- quantreg::rq(y.star ~ xp.star, tau = 0.5)

      # Bootstrap future values (use original last p values)
      x.fut[n:(n - p + 1)] <- series[n:(n - p + 1)]
      for (t in (n + 1):(n + h)) {
        x.fut[t] <- m.lad.star$coef[1] + m.lad.star$coef[-1] %*% x.fut[(t - 1):(t - p)] + res.lad[t]
      }
    }

    p.prr[, b] <- x.fut[(n + 1):(n + h)]
  }

  # --------------------
  # Prediction interval
  # --------------------
  lpi.prr <- apply(p.prr, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  upi.prr <- apply(p.prr, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  l.ppr <- upi.prr - lpi.prr

  list(bfor = p.prr, lpi = lpi.prr, upi = upi.prr, len = l.ppr)
}



