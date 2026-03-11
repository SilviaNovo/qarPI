#' Prediction intervals based on PP algorithm
#'
#' Computes a bootstrap predictive-root-based prediction interval using the forward algorithm
#' of Pan and Politis (2016) for an AR(p) model over a forecast horizon \code{h}.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param B Number of bootstrap replicates. Default is 1000.
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
#' This function implements the PP algorithm of Pan and Politis (2016),
#' referred to in their article as Fp (forward bootstrap with predictive residuals).
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
#' out <- pi_PP(series,alpha=0.10)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' @references
#'Pan, L. and Politis, D. N. (2016). Bootstrap prediction intervals for linear, nonlinear and nonpara-
#'metric autoregressions. \emph{Journal of Statistical Planning and Inference}, 177:1–27.
#'
#' @export
pi_PP <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05) {

  # This function is based on D. Politis code, available on his website
  # https://mathweb.ucsd.edu/~politis/DPpublication.html

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
  if (n - p < 2L)
    stop("Not enough observations to compute predictive residuals.", call. = FALSE)

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))

  y   <- series[(p + 1):n]
  x.p <- as.matrix(xp)

  # --------------------
  # Predictive residuals (leave-one-out)
  # --------------------
  resid <- numeric(n - p)

  for (i in 1:(n - p)) {
    xtemp <- x.p[-i, , drop = FALSE]
    ytemp <- y[-i]
    tempdata <- data.frame(ytemp = ytemp, xtemp = xtemp)
    colnames(tempdata) <- c("ytemp", paste0("xtemp.", 1:p))

    tempfit <- stats::lm(ytemp ~ ., data = tempdata)

    newdata <- as.data.frame(as.list(x.p[i, ]))
    colnames(newdata) <- paste0("xtemp.", 1:p)

    resid[i] <- y[i] - stats::predict(tempfit, newdata = newdata)
  }
  res.p <- resid - base::mean(resid)

  # --------------------
  # OLS estimation
  # --------------------
  m.ols <- stats::lm(y ~ xp)

  # --------------------
  # Resampling (forward bootstrap + predictive root)
  # --------------------
  m <- 100  # burn-in
  b.root <- matrix(0, ncol = B, nrow = h)

  for (b in 1:B) {
    set.seed(b)

    u.star <- numeric(n + m)

    # random p-tuple
    position <- sample(1:(n - p + 1), 1)
    u.star[1:p] <- series[position:(position + p - 1)]

    # resample predictive residuals
    res.ols.star <- sample(res.p, n + m, replace = TRUE)

    # forward bootstrap resampling
    for (i in (p + 1):(n + m)) {
      u.star[i] <- m.ols$coef[1] + m.ols$coef[-1] %*% u.star[(i - 1):(i - p)] + res.ols.star[i]
    }

    x.star <- u.star[(m + 1):(n + m)]

    # bootstrap version of coefficients
    if (p == 1) xp.star <- sapply(1:(n - p), function(x) x.star[(x + p - 1):x])
    if (p > 1)  xp.star <- t(sapply(1:(n - p), function(x) x.star[(x + p - 1):x]))
    y.star <- x.star[(p + 1):n]
    m.ols.star <- stats::lm(y.star ~ xp.star)

    # bootstrap predicted values x^* (predicted under bootstrapped coefficients)
    x.star.hat.fut <- numeric(p + h)
    x.star.hat.fut[1:p] <- series[(n - p + 1):n]
    for (t in 1:h) {
      x.star.hat.fut[t + p] <- m.ols.star$coef[1] +
        m.ols.star$coef[-1] %*% x.star.hat.fut[(p + t - 1):t]
    }
    x.star.hat.f <- x.star.hat.fut[(p + 1):(p + h)]

    # bootstrap future values x* (data-generating with predictive residuals)
    res.ols.p <- sample(res.p, h + p, replace = TRUE)
    x.star.fut <- numeric(p + h)
    x.star.fut[1:p] <- series[(n - p + 1):n]
    for (t in 1:h) {
      x.star.fut[t + p] <- m.ols$coef[1] +
        m.ols$coef[-1] %*% x.star.fut[(p + t - 1):t] +
        res.ols.p[t + p]
    }
    x.star.f <- x.star.fut[(p + 1):(p + h)]

    # predictive root
    b.root[, b] <- x.star.f - x.star.hat.f
  }

  # --------------------
  # Prediction interval
  # --------------------
  x.hat.fut <- numeric(p + h)
  x.hat.fut[1:p] <- series[(n - p + 1):n]
  for (t in 1:h) {
    x.hat.fut[t + p] <- m.ols$coef[1] + m.ols$coef[-1] %*% x.hat.fut[(p + t - 1):t]
  }
  x.hat.f <- x.hat.fut[(p + 1):(p + h)]

  q.alpha  <- apply(b.root, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  q.1alpha <- apply(b.root, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  lpi.fp <- x.hat.f + q.alpha
  upi.fp <- x.hat.f + q.1alpha
  l.fp   <- upi.fp - lpi.fp

  list(pfor = x.hat.f, lpi = lpi.fp, upi = upi.fp, len = l.fp)
}
