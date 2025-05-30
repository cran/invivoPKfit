% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_methods_get_fit.R
\name{get_fit.pk}
\alias{get_fit.pk}
\title{Get fits from a `pk` object}
\usage{
\method{get_fit}{pk}(obj, model = NULL, ...)
}
\arguments{
\item{obj}{A [pk] object.}

\item{model}{Optional: Specify one or more of the fitted models for which to
make predictions. If NULL (the default), predictions will be returned for
all of the models in `obj$stat_model`.}

\item{...}{Additional arguments. Not in use.}
}
\value{
A named list of objects of class `optimx`, named for the models in
  `model`. As described in [optimx::optimx()]  If only one model is
  specified, the return value will still be a list, but with only one
  element.
}
\description{
Get the [optimx::optimx()] output from a fitted `pk` object
}
\details{
This function returns the object(s) returned by [optimx::optimx()] for the
specified model(s) and method(s), for a fitted `pk` object. See
[optimx::optimx()] for details. Briefly, an `optimx` object is a `data.frame`
with one row for each method used, and variables that give the optimized
values for each parameter, along with several diagnostic variables (e.g. the
objective function value at the optimized parameter values; the number of
function evaluations/iterations; an integer code describing convergence
status). The object will have attributes `details` (providing
any messages returned by the methods) and `npar` (the number of parameters
optimized).
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
\code{\link{get_hessian.pk}()},
\code{\link{get_tkstats.pk}()},
\code{\link{logLik.pk}()},
\code{\link{predict.pk}()},
\code{\link{residuals.pk}()},
\code{\link{rmse.pk}()},
\code{\link{rsq.pk}()}
}
\author{
Caroline Ring, Gilberto Padilla Mercado
}
\concept{methods for fitted pk objects}
