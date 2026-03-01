#' Prediction intervals based on AR-perc algorithm
#'
#' Computes a bootstrap percentile-based prediction interval using the AR-perc algorithm for
#' AR(p) model over a forecast horizon \code{h}.
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
#'   \item{bfor}{Numeric matrix of bootstrap forecasts with dimension h x B.}
#'   \item{lpi}{Numeric vector of lower bounds (length \code{h}).}
#'   \item{upi}{Numeric vector of upper bounds (length \code{h}).}
#'   \item{len}{Numeric vector of interval lengths (length \code{h}).}
#' }
#'
#' @details
#'
#' This function implements the AR-perc algorithm described in
#' Novo and Sánchez-Sellero (2025).
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
#' out <- pi_AR_perc(series)
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
#' @export
pi_AR_perc <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05, tau = 0.5) {

  # --------------------
  # Argument checks
  # --------------------
  if (missing(series)) stop("`series` is required.", call. = FALSE)
  if (!is.numeric(series) && !stats::is.ts(series)) {
    stop("`series` must be a numeric vector or a `ts` object.", call. = FALSE)
  }
  series <- as.numeric(series)
  if (anyNA(series)) stop("`series` contains NA values.", call. = FALSE)

  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required. Install it with install.packages('quantreg').", call. = FALSE)
  }

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

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))
  y <- series[(p + 1):n]

  m.quantile <- quantreg::rq(y ~ xp, tau = tau)

  # --------------------
  # Resampling
  # --------------------
  p.qar <- matrix(0, ncol = B, nrow = h)

  for (b in 1:B) {
    set.seed(b)

    # Resample residuals
    r.qar <- sample(m.quantile$residuals, h, replace = TRUE)

    # Multiplier bootstrap for coefficients
    w <- stats::rexp(n - p)
    m.quantile.star <- quantreg::rq(y ~ xp, tau = tau, weights = w)

    # Bootstrap future values
    x.fut <- numeric(p + h)
    x.fut[1:p] <- series[(n - p + 1):n]

    for (t in 1:h) {
      x.fut[t + p] <- m.quantile.star$coef[1] +
        m.quantile.star$coef[-1] %*% x.fut[(p + t - 1):t] +
        r.qar[t]
    }

    p.qar[, b] <- x.fut[(p + 1):(h + p)]
  }

  # --------------------
  # Prediction interval
  # --------------------
  lpi.qar <- apply(p.qar, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  upi.qar <- apply(p.qar, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  l.qar <- upi.qar - lpi.qar

  list(bfor = p.qar, lpi = lpi.qar, upi = upi.qar, len = l.qar)
}


