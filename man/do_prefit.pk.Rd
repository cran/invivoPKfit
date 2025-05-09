% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_methods_do_prefit.R
\name{do_prefit.pk}
\alias{do_prefit.pk}
\title{Do pre-fitting}
\usage{
\method{do_prefit}{pk}(obj, ...)
}
\arguments{
\item{obj}{A `pk` object}

\item{...}{Additional arguments. Not in use.}
}
\value{
The same `pk` object, but with a new element `prefit`, containing the
  results of pre-fit calculations and checks for each model and for the error
  model.
}
\description{
Do pre-fit calculations and checks
}
\details{
This function does the following:

- Based on the error model in `stat_error_model` and the pre-processed data, determines the number of residual standard deviations ("sigmas") hyperparameters to be estimated.
- Determines which "sigma" hyperparameter corresponds to each observation in the data.
- Calculates lower/upper bounds and starting guesses for each "sigma" hyperparameter
- For each model in `stat_model`, calls its `params_fun`, the function that, based on the data, determines whether to optimize each model parameter, and calculates lower/upper bounds and starting guesses for each model parameter to be optimized. Only non-excluded observations are passed to each model's `params_fun`.


Lower bounds for each "sigma" hyperparameter are set to `sqrt(.Machine$double_eps)`.

Upper bounds for each "sigma" hyperparameter are calculated as the standard
deviation of observations in the corresponding error SD group (see
[combined_sd()]), with any specified transformations applied
(dose-normalization and/or log10-transformation). If the combined SD is
non-finite or less than the sigma lower bound, then the maximum concentration
is used as an upper bound; if this still returns a non-finite value or a
value less than the lower bound, then a constant value of 1000 is
substituted.

The starting guess for each "sigma" hyperparameter is one-tenth of the upper bound.
}
\author{
Caroline Ring
}
