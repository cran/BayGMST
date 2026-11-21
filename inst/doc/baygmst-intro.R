## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)
has_cmdstan <- requireNamespace("instantiate", quietly = TRUE) &&
  instantiate::stan_cmdstan_exists()

## ----setup--------------------------------------------------------------------
library(BayGMST)

## ----load-data----------------------------------------------------------------
proxy_matrix <- read.csv(
  system.file("extdata", "example_proxy_matrix.csv", package = "BayGMST")
)
forcing_raw <- read.csv(
  system.file("extdata", "example_forcing.csv", package = "BayGMST")
)
instrumental_raw <- read.csv(
  system.file("extdata", "example_instrumental_temperature.csv", package = "BayGMST")
)
# read.csv() mangles "Time"/"Anomaly (deg C)" headers via make.names(); give
# them stable names instead of relying on the mangled spelling.
colnames(instrumental_raw)[1:2] <- c("Time", "Anomaly")

years <- forcing_raw$year
head(proxy_matrix)

## ----reduce-proxies-----------------------------------------------------------
calib_years <- instrumental_raw$Time

rp <- reduce_proxies(
  proxy_matrix = proxy_matrix[, -1],
  years        = proxy_matrix$year,
  temp_calib   = instrumental_raw$Anomaly,
  calib_years  = calib_years,
  method       = "PCR",
  # `chunk` must be large enough that every segment (including the last,
  # shortest one) still fully contains the calibration window (1850-2000,
  # 151 years). This toy example only has 301 years total, so we use a
  # single segment (chunk >= length(years)) rather than the package
  # default of 250, which would create a too-short trailing segment here.
  chunk        = nrow(proxy_matrix)
)
rp

## ----transform-forcings-------------------------------------------------------
forcing <- transform_forcings(
  G = forcing_raw$CO2,
  V = forcing_raw$volcanic,
  S = forcing_raw$solar
)

## ----align--------------------------------------------------------------------
instrumental_T <- instrumental_raw$Anomaly[match(years, instrumental_raw$Time)]

## ----fit, eval = has_cmdstan--------------------------------------------------
fit <- fit_baygmst(
  proxy          = rp,
  instrumental_T = instrumental_T,
  forcing_G      = forcing$G,
  forcing_V      = forcing$V,
  forcing_S      = forcing$S,
  years          = years,
  chains         = 2,
  parallel_chains = 2,
  iter_warmup    = 300,
  iter_sampling  = 1000,
  seed           = 1
)
summary(fit)

## ----plot-fit, eval = has_cmdstan && requireNamespace("patchwork", quietly = TRUE), fig.height = 6----
plot(fit)

## ----plot-reconstruction, eval = has_cmdstan----------------------------------
plot_reconstruction(fit)

## ----plot-trace, eval = has_cmdstan-------------------------------------------
plot_trace(fit)

## ----plot-posteriors, eval = has_cmdstan--------------------------------------
plot_posterior_densities(fit)

## ----cv, eval = has_cmdstan---------------------------------------------------
cv <- cv_baygmst(
  proxy          = rp,
  instrumental_T = instrumental_T,
  forcing_G      = forcing$G,
  forcing_V      = forcing$V,
  forcing_S      = forcing$S,
  years          = years,
  nfold          = 3,
  chains         = 2,
  parallel_chains = 2,
  iter_warmup    = 300,
  iter_sampling  = 1000,
  seed           = 1
)
cv
plot_cv(cv)

## ----byo-composite------------------------------------------------------------
set.seed(1)
composite_ens <- list(
  ages      = seq(250, -50, by = -1),               # 1700-2000 CE, in yr BP
  composite = matrix(rnorm(301 * 50, sd = 0.2), 301, 50) +
    seq(-0.2, 0.4, length.out = 301)                # 50 ensemble members
)

proxy_input <- as_baygmst_proxy(composite_ens)
proxy_input

