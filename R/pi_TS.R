#' Prediction intervals based on the TS algorithm
#'
#' Computes bootstrap percentile-based prediction intervals using the backward
#' algorithm of Thombs and Schucany (1990) for an AR(p) model over a
#' forecast horizon \code{h}.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param B Number of bootstrap replicates. Default is 1000.
#' @param alpha Significance level. 1-\code{alpha} is the nominal coverage level. Default is 0.05.
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
#' This function implements the backward bootstrap algorithm proposed by
#' Thombs and Schucany (1990).
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
#' out <- pi_TS(series,h=4)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' @references
#' Thombs, L. A. and Schucany, W. R. (1990). Bootstrap prediction intervals for autoregression.
#' \emph{Journal of the American Statistical Association}, 85(410):486–492.
#'
#' @export

pi_TS <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05) {

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

  if (!is.numeric(p) || length(p) != 1L || is.na(p))
    stop("`p` must be a single number.", call. = FALSE)
  if (p != as.integer(p) || p < 1L)
    stop("`p` must be a positive integer.", call. = FALSE)

  if (!is.numeric(h) || length(h) != 1L || is.na(h))
    stop("`h` must be a single number.", call. = FALSE)
  if (h != as.integer(h) || h < 1L)
    stop("`h` must be a positive integer.", call. = FALSE)

  if (!is.numeric(B) || length(B) != 1L || is.na(B))
    stop("`B` must be a single number.", call. = FALSE)
  if (B != as.integer(B) || B < 2L)
    stop("`B` must be a positive integer >= 2.", call. = FALSE)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha))
    stop("`alpha` must be a single number.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).", call. = FALSE)

  if (n <= p)
    stop("`series` length must be greater than `p`.", call. = FALSE)
  if (n - 2L * p <= 0L)
    stop("Need `n - 2*p > 0` for the residual rescaling used in this function.", call. = FALSE)

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))
  y <- series[(p + 1):n]

  m.ols <- stats::lm(y ~ xp)

  # 1) Backward residuals
  r.back <- numeric(n - p)
  for (i in 1:(n - p)) {
    r.back[i] <- series[i] - m.ols$coef[1] - m.ols$coef[-1] %*% series[(i + 1):(i + p)]
  }
  r.back.res <- sqrt((n - p) / (n - 2 * p)) * (r.back - base::mean(r.back))

  # 2) Forward residuals
  r.ols <- m.ols$residuals
  r.ols.res <- sqrt((n - p) / (n - 2 * p)) * (r.ols - base::mean(r.ols))

  # --------------------
  # Bootstrap
  # --------------------
  p.ts <- matrix(0, ncol = B, nrow = h)
  for (b in 1:B) {
    set.seed(b)

    # resample backward residuals
    res.back <- sample(r.back.res, n - p, replace = TRUE)

    # backward bootstrap sample
    x.fut <- numeric(n + h)
    x.fut[n:(n - p + 1)] <- series[n:(n - p + 1)]
    for (i in (n - p):1) {
      x.fut[i] <- m.ols$coef[1] + m.ols$coef[-1] %*% x.fut[(i + 1):(i + p)] + res.back[i]
    }

    # refit on bootstrap sample
    if (p == 1) xp.star <- sapply(1:(n - p), function(x) x.fut[(x + p - 1):x])
    if (p > 1)  xp.star <- t(sapply(1:(n - p), function(x) x.fut[(x + p - 1):x]))
    y.star <- x.fut[(p + 1):n]
    m.ols.star <- stats::lm(y.star ~ xp.star)

    # resample forward residuals
    res.for <- sample(r.ols.res, n + h, replace = TRUE)

    # forward bootstrap for future
    for (i in (n + 1):(n + h)) {
      x.fut[i] <- m.ols.star$coef[1] + m.ols.star$coef[-1] %*% x.fut[(i - 1):(i - p)] + res.for[i]
    }
    p.ts[, b] <- x.fut[(n + 1):(n + h)]
  }

  # --------------------
  # Prediction interval
  # --------------------
  lpi.ts <- apply(p.ts, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  upi.ts <- apply(p.ts, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  len.ts <- upi.ts - lpi.ts

  list(bfor = p.ts, lpi = lpi.ts, upi = upi.ts, len = len.ts)
}
