#' Prediction intervals based on the CB algorithm
#'
#' Computes bootstrap percentile-based prediction intervals using the
#' conditional bootstrap algorithm of Cao et al. (1997) for
#' an AR(p) model over a forecast horizon \code{h}.
#'
#' @param series Numeric vector or ts object with time series values.
#' @param p Positive integer indicating the autoregressive order. Default is 1.
#' @param h Positive integer indicating the prediction horizon. Default is 3.
#' @param B Number of bootstrap replicates. Default is 1000.
#' @param alpha Significance level. 1-\code{alpha} is the nominal coverage level. Default is 0.05.
#' @param method Estimation method. One among "OLS" and "LAD". Default is "OLS".
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
#' This function implements the conditional bootstrap algorithm described in
#' Cao et al. (1997).
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
#' out <- pi_CB(series,h=4,alpha=0.01)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' @references
#' Cao, R., Febrero-Bande, M., González-Manteiga, W., Prada-Sánchez, J., and García-Jurado, I. (1997).
#' Saving computer time in constructing consistent bootstrap prediction intervals for autoregressive
#' processes. \emph{Communications in Statistics-Simulation and Computation}, 26(3):961–978.
#'
#' @export
pi_CB <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05,
                  method = c("OLS", "LAD")) {

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
  if (n - 2L * p <= 0L) {
    stop("Need `n - 2*p > 0` for the residual rescaling used in this function.", call. = FALSE)
  }

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
  p.cb <- matrix(0, ncol = B, nrow = h)

  for (b in 1:B) {
    set.seed(b)

    if (method == "OLS") {
      res.ols <- sample(r.ols.res, h, replace = TRUE)

      x.fut <- numeric(p + h)
      x.fut[1:p] <- series[(n - p + 1):n]
      for (t in 1:h) {
        x.fut[t + p] <- m.ols$coef[1] + m.ols$coef[-1] %*% x.fut[(p + t - 1):t] + res.ols[t]
      }
    }

    if (method == "LAD") {
      res.lad <- sample(r.lad, h, replace = TRUE)

      x.fut <- numeric(p + h)
      x.fut[1:p] <- series[(n - p + 1):n]
      for (t in 1:h) {
        x.fut[t + p] <- m.lad$coef[1] + m.lad$coef[-1] %*% x.fut[(p + t - 1):t] + res.lad[t]
      }
    }

    p.cb[, b] <- x.fut[(p + 1):(h + p)]
  }

  # --------------------
  # Prediction interval
  # --------------------
  lpi.cb <- apply(p.cb, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  upi.cb <- apply(p.cb, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  l.cb <- upi.cb - lpi.cb

  list(bfor = p.cb, lpi = lpi.cb, upi = upi.cb, len = l.cb)
}


