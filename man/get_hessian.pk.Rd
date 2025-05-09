% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_methods_get_hessian.R
\name{get_hessian.pk}
\alias{get_hessian.pk}
\title{Get Hessian matrixes}
\usage{
\method{get_hessian}{pk}(obj, model = NULL, method = NULL, suppress.messages = TRUE, ...)
}
\arguments{
\item{obj}{A [pk] object}

\item{model}{Optional: Specify one or more of the fitted models whose
coefficients to return. If NULL (the default), coefficients will be returned
for all of the models in `obj$stat_model`.}

\item{method}{Optional: Specify one or more of the [optimx::optimx()] methods
whose coefficients to return. If NULL (the default), coefficients will be
returned for all of the models in `obj$settings_optimx$method`.}

\item{suppress.messages}{Logical. `TRUE` (the default) to suppress informative
messages. `FALSE` to see them.}

\item{...}{Additional arguments. Not in use right now.}
}
\value{
A dataframe with one row for each `data_group`, `model` and `method`.
 The remaining column is a `list` column containing the Hessian for each row.
}
\description{
Extract Hessian matrixes from a fitted `pk` object
}
\details{
This function computes a numerical approximation to the model Hessian for each
data group and each model in a fitted `pk` object. The Hessian is the matrix
of second derivatives of the model objective function with respect to each
model parameter. Here, the objective function is the negative log-likelihood
implemented in [log_likelihood()], evaluated jointly across the data that was
used to fit the model.
}
\references{
Gill J, King G. (2004) What to Do When Your Hessian is Not
 Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
 Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
}
\seealso{
Other methods for fitted pk objects: 
\code{\link{AAFE.pk}()},
\code{\link{AFE.pk}()},
\code{\link{AIC.pk}()},
\code{\link{BIC.pk}()},
\code{\link{coef.pk}()},
\code{\link{coef_sd.pk}()},
\code{\link{eval_tkstats.pk}()},
\code{\link{get_fit.pk}()},
\code{\link{get_tkstats.pk}()},
\code{\link{logLik.pk}()},
\code{\link{predict.pk}()},
\code{\link{residuals.pk}()},
\code{\link{rmse.pk}()},
\code{\link{rsq.pk}()}
}
\author{
Caroline Ring and Gilberto Padilla Mercado
}
\concept{methods for fitted pk objects}
