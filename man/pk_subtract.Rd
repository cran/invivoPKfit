% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_methods_add_pkproto.R
\name{pk_subtract}
\alias{pk_subtract}
\title{Subtract a `pkproto` object from a `pk` object}
\usage{
pk_subtract(pkproto_obj, pk_obj, objectname)
}
\arguments{
\item{pkproto_obj}{The `pkproto` object to be subtracted}

\item{pk_obj}{The `pk` object to which the `pkproto` object is to be subtracted}

\item{objectname}{The object name}
}
\value{
The `pk` object modified by the subtraction.
}
\description{
This is the S3 generic method.
}
\seealso{
[pk_subtract.pk_stat_model()] for the method for subtracting `pk_stat_model` objects
  (from `stat_model()`)
}
