% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_starts_1comp.R
\name{get_starts_1comp}
\alias{get_starts_1comp}
\title{Get starting values for 1-compartment model}
\usage{
get_starts_1comp(data, par_DF)
}
\arguments{
\item{data}{The data set to be fitted (e.g. the result of [preprocess_data()])}

\item{par_DF}{A `data.frame` with the following variables (e.g., as produced by [get_params_1comp()])
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
Derive starting values for 1-compartment model parameters from available data
}
\details{
This function is called internally by [get_params_1comp()] and should
generally not be called directly by the user.

The full set of model parameters for the 1-compartment model includes `Vdist`,
`kelim`, `kgutabs`, `Fgutabs`, and `Rblood2plasma`. Whether each one can be
estimated from the data depends on what routes of administration are included
in the data.

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

# Starting value for `kelim`

If IV data exist, then only IV data are used to derive starting estimates for
`kelim`, even if oral data also exist.

If only oral data exist, then the oral data are used to derive a starting
estimate for `kelim`.

Whichever data set is used (IV or oral), the starting value for `kelim` is
derived by assuming that the range of observed time values in the data set
spans two elimination half-lives. This implies that the elimination half-life
is equal to the midpoint of observed time values, and that the starting value
for the elimination time constant `kelim` is therefore `log(2)` divided by the
midpoint of observed time values.

Of course, this assumption is unlikely to be correct. However, we hope that it
will yield a starting guess for `kelim` that is at least on the right order of
magnitude.

# Starting value for `Vdist`

If IV data exist, then only IV data are used to derive a starting estimate
for `Vdist`.

This starting estimate is derived by assuming that the IV data obey a
one-compartment model, which means that when concentrations are
dose-normalized and log10-transformed and plotted against time, they will
follow a straight line with slope `-kelim`.

First, concentrations are dose-normalized by dividing them by their corresponding
doses. Then the normalized concentrations are log10-transformed.

From all observations at the earliest observed time point in the data set
(call it `tmin`), the median of the dose-normalized, log10-transformed
concentrations is calculated; call it `C_tmin`. (The median is used, rather
than the mean, in an attempt to be more robust to outliers.)

If the earliest observed time point is not at time = 0, then the
dose-normalized, log10-transformed concentration at time = 0 is extrapolated
by drawing a straight line with slope `-kelim` back from `C_tmin`, where the
value of `kelim` is the starting value derived as in the previous section.

This extrapolated concentration at time t = 0 is called `A_log10`. `A_log10`
represents the expected body concentration immediately after IV injection of
a unit dose (under the assumption that TK obeys a one-compartment model).

Then, the volume of distribution `Vdist` is derived as `1/(10^A_log10)`. In
other words, `Vdist` is the volume that would be required to produce a
concentration equal to `A_log10` after injecting a unit dose.

(No starting value for `Vdist` can be derived with only oral data,
but none is needed, because with only oral data, `Vdist` will not be estimated
from the data).

# Starting value for `kgutabs`

If oral data exist (whether or not IV data also exist), then the oral data
are used to derive a starting value for `kgutabs`.

First, concentrations are dose-normalized by dividing them by their corresponding
doses. Then the normalized concentrations are log10-transformed.

The time of peak concentration (`tmax`), and the median (normalized,
log-transformed) peak concentration (`Cmax_log10`), are identified using [get_peak()].

As a very rough guess,`tmax` is assumed to
occur at one absorption half-life. Under this assumption, `kgutabs` is equal
to `log(2)/tmax`, and this is taken as the starting value.

# Starting value for `Fgutabs_Vdist`

If any oral data exist (whether or not IV data also exist), then the oral data
are used to derive a starting value for `Fgutabs_Vdist`.

If the kinetics obey a one-compartment model, then if concentrations are
dose-normalized, log-transformed, and plotted vs. time, then at late time
points (after concentration has peaked), the concentration vs. time
relationship will approach a straight line with slope `-kelim`.

If this straight line is extrapolated back to time 0, then the resulting
intercept (call it `A`), expressed on the natural scale, is equal to
`Fgutabs_Vdist * kgutabs/(kgutabs-kelim)`. See
https://www.boomer.org/c/p4/c09/c0902.php .

Roughly, we approximate `A` on the log10 scale by extrapolating back from the peak along a
straight line with slope `-kelim`, using the previously-derived starting
value for `kelim`. So `log10(A) = Cmax_log10 + kelim*tmax`.

Using the previously-derived starting values for `kgutabs` and `kelim`, then,
the starting value for `Fgutabs_Vdist` can be derived as `A * (kgutabs-kelim)/kgutabs`.

# Starting value for `Fgutabs`

If both oral and IV data exist, then the derived starting values for `Vdist`
(from the IV data) and `Fgutabs_Vdist` (from the oral data) are multiplied to
yield a derived starting value for `Fgutabs`.

#Starting value for `Rblood2plasma`

The starting value for `Rblood2plasma` is always set at a constant 1.
}
\seealso{
Other 1-compartment model functions: 
\code{\link{auc_1comp}()},
\code{\link{auc_1comp_cl}()},
\code{\link{cp_1comp}()},
\code{\link{cp_1comp_cl}()},
\code{\link{get_params_1comp}()},
\code{\link{get_params_1comp_cl}()},
\code{\link{get_params_1comp_fup}()},
\code{\link{get_starts_1comp_cl}()},
\code{\link{get_starts_1comp_fup}()}

Other get_starts functions: 
\code{\link{get_starts_1comp_cl}()},
\code{\link{get_starts_1comp_fup}()},
\code{\link{get_starts_2comp}()},
\code{\link{get_starts_flat}()}

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
\code{\link{get_starts_1comp_cl}()},
\code{\link{get_starts_1comp_fup}()},
\code{\link{get_starts_2comp}()},
\code{\link{get_starts_flat}()},
\code{\link{tkstats_2comp}()},
\code{\link{transformed_params_2comp}()}
}
\author{
Caroline Ring
}
\concept{1-compartment model functions}
\concept{built-in model functions}
\concept{get_starts functions}
