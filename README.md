
<!-- README.md is generated from README.Rmd. Please edit README.Rmd -->

# qarPI

<!-- badges: start -->
<!-- Si activas GitHub Actions con usethis::use_github_action_check_standard(),
     descomenta el badge de R-CMD-check y ajusta USUARIO/REPO -->

<!--[![R-CMD-check](https://github.com/SilviaNovo/qarPI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SilviaNovo/qarPI/actions/workflows/R-CMD-check.yaml) -->
[![License: GPL
v3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2512.22018-b31b1b.svg)](https://doi.org/10.48550/arXiv.2512.22018)
<!-- badges: end -->

The package provides prediction intervals (PIs) under two model classes:

- Classical autoregressive models: **AR(p)**
- Quantile autoregressive models: **QAR(p)**

A central feature of the package is the implementation of the new
bootstrap algorithms **AR-perc**, **AR-proot**, **QAR-perc**, and
**QAR-proot**, introduced in the article *“Prediction Intervals for
Quantile Autoregression”*
([arXiv:2512.22018](https://doi.org/10.48550/arXiv.2512.22018)).

In addition, the package includes several established bootstrap and
analytical procedures (Box-Jenkins) from the existing literature,
allowing users to replicate and compare alternative PI constructions
within a unified framework.

Within each model class, bootstrap PIs may therefore be obtained via:

- Percentile-based bootstrap procedures
- Predictive-root-based bootstrap procedures

------------------------------------------------------------------------

### Installation

``` r
# install the package
devtools::install_github("SilviaNovo/qarPI")


#load the package
library(qarPI)
```

------------------------------------------------------------------------

### 1. PI methods for AR(p) models

#### Percentile-based bootstrap PIs

Bootstrap procedures constructed under the AR(p) framework using the
percentile principle

- **AR-perc** — Bootstrap algorithm based on quantile estimation and
  bootstrap multipliers (Novo & Sanchez-Sellero, 2025) → `pi_AR_perc()`

``` r
set.seed(123)
y <- arima.sim(model = list(ar = 0.6), n = 300)

# AR-perc PI
pi <- pi_AR_perc(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.849219 1.828965 3.678184
    ## 2 -2.122772 2.202761 4.325532
    ## 3 -2.171773 2.409590 4.581363

- **TS** — Backward bootstrap (Thombs & Schucany, 1990)  
  → `pi_TS()`

``` r
# Thombs & Schucany PI
pi <- pi_TS(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.809166 1.961925 3.771091
    ## 2 -2.062066 2.340589 4.402655
    ## 3 -2.116774 2.275783 4.392557

- **CB** — Conditional bootstrap (Cao et al., 1997)  
  → `pi_CB()`

``` r
# Cao et al. PI
pi <- pi_CB(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.819219 1.889094 3.708313
    ## 2 -2.085047 2.205806 4.290853
    ## 3 -2.068545 2.460602 4.529146

- **PRR** — Forward bootstrap (Pascual, Romo & Ruiz, 2004) using
  ordinary least squares (OLS) or least absolute deviations (LAD)
  estimation → `pi_PRR()`

``` r
# Pascual et al. PI. OLS
pi <- pi_PRR(y, p = 1, h = 3, method="OLS")

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.866294 1.906787 3.773081
    ## 2 -2.003852 2.309888 4.313740
    ## 3 -2.116450 2.481027 4.597478

``` r
# Pascual et al. PI. LAD
pi <- pi_PRR(y, p = 1, h = 3, method="LAD")

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.883711 1.882117 3.765828
    ## 2 -1.993665 2.283575 4.277240
    ## 3 -2.159428 2.409537 4.568965

------------------------------------------------------------------------

#### Predictive-root-based bootstrap PIs

Bootstrap intervals constructed using the predictive root approach

- **AR-proot** — Bootstrap algorithm based on quantile estimation,
  bootstrap multipliers and predictive residuals (Novo &
  Sanchez-Sellero, 2025) → `pi_AR_proot()`

``` r
# AR-proot PI
pi <- pi_AR_proot(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.778688 1.944880 3.723568
    ## 2 -2.063616 2.239408 4.303024
    ## 3 -2.043930 2.589753 4.633683

- **PP** — Forward bootstrap with predictive residuals (Pan & Politis,
  2016, denoted by the authors as *Fp*). → `pi_PP()`

``` r
# Pan & Politis PI
pi <- pi_PP(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.848580 1.952997 3.801577
    ## 2 -2.159115 2.209435 4.368550
    ## 3 -2.162428 2.142813 4.305242

#### Analytical PIs: Box-Jenkins

- **BJ** — Box–Jenkins type analytical PIs assuming Gaussian
  innovations  
  → `pi_BJ()`

``` r
# Classical Box–Jenkins PI
pi <- pi_BJ(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -2.023091 1.737312 3.760403
    ## 2 -2.205069 2.090839 4.295908
    ## 3 -2.232970 2.213495 4.446464

------------------------------------------------------------------------

------------------------------------------------------------------------

### 2. PI methods for QAR(p) models

Prediction intervals constructed under the QAR(p) framework

------------------------------------------------------------------------

#### Percentile-based bootstrap PIs

- **QAR-perc** — Bootstrap algorithm based on quantile estimation and
  bootstrap multipliers (Novo & Sanchez-Sellero, 2025) → `pi_QAR_perc()`

``` r
# QAR-perc PI
pi <- pi_QAR_perc(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.858576 1.860805 3.719380
    ## 2 -2.137675 2.037932 4.175606
    ## 3 -2.113105 2.394030 4.507136

- **X** — Bootstrap algorithm suggested by Xiao (2012)  
  → `pi_X()`

``` r
# Xiao PI
pi <- pi_X(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.836173 1.916079 3.752252
    ## 2 -2.066253 2.279619 4.345872
    ## 3 -2.163826 2.320773 4.484599

------------------------------------------------------------------------

#### Predictive-root-based bootstrap PIs

- **QAR-proot** — Bootstrap algorithm based on quantile estimation,
  bootstrap multipliers and predictive residuals (Novo &
  Sanchez-Sellero, 2025) → `pi_QAR_proot()`

``` r
# QAR-proot PI
pi <- pi_QAR_proot(y, p = 1, h = 3)

t_PI <- data.frame(
  Lower = pi$lpi,
  Upper = pi$upi,
  Length = pi$len
)
t_PI
```

    ##       Lower    Upper   Length
    ## 1 -1.826982 1.894856 3.721838
    ## 2 -2.051728 2.082905 4.134633
    ## 3 -2.047210 2.464404 4.511614

------------------------------------------------------------------------

## Summary Structure

``` text
qarPI
│
├── AR(p)
│   ├── Percentile-based Bootstrap
│   │   ├── AR-perc → pi_AR_perc()
│   │   ├── TS → pi_TS()
│   │   ├── CB → pi_CB()
│   │   └── PRR → pi_PRR()
│   │
│   ├── Predictive-root-based Bootstrap
│   │   ├── AR-proot → pi_AR_proot()
│   │   └── PP (Fp) → pi_PP()
│   │
│   └── Analytical
│       └── BJ → pi_BJ()
│
└── QAR(p)
    ├── Percentile-based Bootstrap
    │   ├── QAR-perc → pi_QAR_perc()
    │   └── X → pi_X()
    │
    └── Predictive-root-based Bootstrap
        └── QAR-proot → pi_QAR_proot()
```

## References

- Box, G. E. & Jenkins, G. M. (1976). *Time Series Analysis: Forecasting
  and Control*. Holden–Day, San Francisco.

- Cao, R., et al. (1997). Saving computer time in constructing
  consistent bootstrap prediction intervals for autoregressive
  processes. *Communications in Statistics-Simulation and Computation*,
  26(3):961–978. <https://doi.org/10.1080/03610919708813420>

- Novo, S., & Sánchez-Sellero, C. A. Prediction intervals for quantile
  autoregression. <https://doi.org/10.48550/arXiv.2512.22018>

- Pan, L., & Politis, D. N. (2016). Bootstrap prediction intervals for
  linear, nonlinear and nonparametric autoregressions. *Journal of
  Statistical Planning and Inference*.
  <https://doi.org/10.1016/j.jspi.2014.10.003>

- Pascual, L., Romo, J., & Ruiz, E. (2004). Bootstrap predictive
  inference for ARIMA processes. *Journal of Time Series Analysis*.
  25(4):449–465. <https://doi.org/10.1111/j.1467-9892.2004.01713.x>

- Thombs, L. A., & Schucany, W. R. (1990). Bootstrap prediction
  intervals for autoregression. *Journal of the American Statistical
  Association*, 85(410):486–492.
  <https://doi.org/10.1080/01621459.1990.10476225>

- Xiao, Z. (2012). Time series quantile regressions. Handbook of
  Statistics. Time Series Analysis: Methods and Applications,
  30:213–257. <https://doi.org/10.1016/B978-0-444-53858-1.00009-0>
