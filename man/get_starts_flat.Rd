% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_starts_flat.R
\name{get_starts_flat}
\alias{get_starts_flat}
\title{Get starting values for flat model}
\usage{
get_starts_flat(data, par_DF)
}
\arguments{
\item{data}{The data set to be fitted (e.g. the result of [preprocess_data()])}

\item{par_DF}{A `data.frame` with the following variables (e.g., as produced by [get_params_flat()])
- `param_name`: Character: Names of the model parameters
- `param_units`: Character: Units of the model parameters
- `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
- `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
-`lower_bounds`: Numeric: The lower bounds for each parameter
- `upper_bounds`: Numeric: The upper bounds for each parameter}
}
\value{
The same `data.frame` as `par_DF`, with an additional variable
 `starts` containing the derived starting value for each parameter. If a
 parameter cannot be estimated from the available data, then its starting value
 will be `NA_real_`
}
\description{
Derive starting values for flat model parameters from available data
}
\details{
This function is called internally by [get_params_1comp()] and should
generally not be called directly by the user.

The full set of model parameters for the flat model includes
`Vdist`,`Fgutabs`, and `Rblood2plasma`. Whether each one can be estimated from
the data depends on what routes of administration are included in the data.

The numerical optimizer requires starting guesses for the value of each
parameter to be estimated from the data. Default starting guesses are derived from the available data.

These are intended to be *very* rough starting guesses, so the algorithm here
is extremely naive. This function is not itself intended to produce valid
estimates for any of the model parameters, and it is highly unlikely to do so.

The derivation process is as follows.

First, data are filtered to exclude any non-detects.

Then, data are split by route of administration, into an IV data set and an oral data
set. (It is possible that either IV or oral data may not be
available for a chemical.)

If IV data exist, then only IV data are used to derive a starting estimate for
`Vdist`. Concentrations are dose-normalized (divided by their corresponding
dose) and log10-transformed. The mean dose-normalized, log10-transformed
concentration is calculated (call it `Cmean_log10`). `Vdist` starting value is
then derived as `1/(10^Cmean_log10)` .

If any oral data exist (whether or not IV data also exist), then the oral data
are used to derive a starting value for `Fgutabs_Vdist`. Concentrations are dose-normalized (divided by their corresponding
dose) and log10-transformed. The mean dose-normalized, log10-transformed
concentration is calculated (call it `Cmean_log10`). `Fgutabs_Vdist` starting value is
then set equal to `10^Cmean_log10` .

#Starting value for `Rblood2plasma`

If both blood and plasma data are available, then the starting value for `Rblood2plasma` is derived as follows.

If IV data are available for both blood and plasma, then the starting value
for `Rblood2plasma` is derived as the ratio of `Vdist` for blood data and
`Vdist` for plasma data.

If oral data, but not IV data, are available for both blood and plasma, then
the starting value for `Rblood2plasma` is derived as the ratio of
`Fgutabs_Vdist` for plasma data and `Fgutabs_Vdist` for blood data.

If only blood data or only plasma data are available, then the starting value for `Rblood2plasma` is set at a constant 1.
}
\seealso{
Other flat model functions: 
\code{\link{auc_flat}()},
\code{\link{cp_flat}()},
\code{\link{get_params_flat}()}

Other get_starts functions: 
\code{\link{get_starts_1comp}()},
\code{\link{get_starts_1comp_cl}()},
\code{\link{get_starts_1comp_fup}()},
\code{\link{get_starts_2comp}()}

Other built-in model functions: 
\code{\link{auc_1comp}()},
\code{\link{auc_1comp_cl}()},
\code{\link{auc_2comp}()},
\code{\link{auc_flat}()},
\code{\link{cp_1comp}()},
\code{\link{cp_1comp_cl}()},
\code{\link{cp_2comp}()},
\code{\link{cp_2comp_dt}()},
\code{\link{cp_flat}()},
\code{\link{get_params_1comp}()},
\code{\link{get_params_1comp_cl}()},
\code{\link{get_params_1comp_fup}()},
\code{\link{get_params_2comp}()},
\code{\link{get_params_flat}()},
\code{\link{get_starts_1comp}()},
\code{\link{get_starts_1comp_cl}()},
\code{\link{get_starts_1comp_fup}()},
\code{\link{get_starts_2comp}()},
\code{\link{tkstats_2comp}()},
\code{\link{transformed_params_2comp}()}
}
\author{
Caroline Ring
}
\concept{built-in model functions}
\concept{flat model functions}
\concept{get_starts functions}
