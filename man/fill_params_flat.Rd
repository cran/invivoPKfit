% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fill_params_flat.R
\name{fill_params_flat}
\alias{fill_params_flat}
\title{Fill parameters for flat model}
\usage{
fill_params_flat(params)
}
\arguments{
\item{params}{Named numeric vector of parameters for the flat model}
}
\value{
A named numeric vector of parameters, with any flat model
  parameters not present in `params` filled with `NA_real_`. If any two of
  `Fgutabs`, `Vdist`, and `Fgutabs_Vdist` were present in `params`, the third will
  be imputed to agree with the other two.
}
\description{
Fill parameters for flat model
}
\author{
Caroline Ring
}
