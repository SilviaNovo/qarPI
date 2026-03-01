#' Prediction intervals based on X algorithm
#'
#' Computes a bootstrap percentile-based prediction interval using the algorithm suggested by
#' Xiao (2012) for a QAR(p) model over a forecast horizon \code{h}.
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
#' This function implements the X algorithm suggested by Xiao (2012).
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
#' out <- pi_X(series)
#'
#' out$lpi
#' out$upi
#' out$len
#'
#' # Simulate QAR(1) process
#' n <- 100
#' e<-runif(m)
#' x <- numeric(m)
#' x[c(1,2)] <- rt(2,3)
#'
#' for (t in 3:m) {x[t]<-0.3*x[t-1]+0.7*e[t]*x[t-2]+qt(e[t],3)}
#' series2<-ts(x[301:m])
#'
#' # Compute prediction interval
#'
#' out2<-pi_X(series2,h=4)
#' out2$lpi
#' out2$upi
#' out2$len
#'
#' @references
#'Xiao, Z. (2012). Time series quantile regressions. Handbook of Statistics. Time Series Analysis:
#'Methods and Applications, 30:213–257.
#'
#'
#'@export
pi_X <- function(series, p = 1, h = 3, B = 1000, alpha = 0.05) {

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

  # --------------------
  # Series transformation in response+predictor format
  # --------------------
  if (p == 1) xp <- sapply(1:(n - p), function(x) series[(x + p - 1):x])
  if (p > 1)  xp <- t(sapply(1:(n - p), function(x) series[(x + p - 1):x]))
  y <- series[(p + 1):n]

  # --------------------
  # Resampling
  # --------------------
  p.qar <- matrix(0, ncol = B, nrow = h)
  for (b in 1:B) {
    set.seed(b)

    # Bootstrap future values
    x.fut <- numeric(p + h)
    x.fut[1:p] <- series[(n - p + 1):n]

    for (t in 1:h) {
      q <- stats::runif(1)
      m.quantile.star <- quantreg::rq(y ~ xp, tau = q)
      x.fut[t + p] <- m.quantile.star$coef[1] +
        m.quantile.star$coef[-1] %*% x.fut[(p + t - 1):t]
    }

    p.qar[, b] <- x.fut[(p + 1):(h + p)]
  }

  # --------------------
  # Prediction interval
  # --------------------
  lpi.qar <- apply(p.qar, 1, sort)[max(1L, floor(alpha / 2 * B)), ]
  upi.qar <- apply(p.qar, 1, sort)[min(B, floor((1 - alpha / 2) * B)), ]

  len.qar <- upi.qar - lpi.qar

  list(bfor = p.qar, lpi = lpi.qar, upi = upi.qar, len = len.qar)
}
