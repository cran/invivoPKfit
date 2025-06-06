#' Toxicokinetic statistics for 1-compartment model
#'
#' Calculate predicted toxicokinetic statistics for a 1-compartment model.
#'
#' # Statistics computed
#'
#' ## Total clearance
#'
#' \deqn{\textrm{CL}_{tot} = k_{elim} + V_{dist}}
#'
#' ## Steady-state plasma concentration for long-term daily dose of 1 mg/kg/day
#'
#' The dosing interval \eqn{\tau = \frac{1}{\textrm{day}}} will be converted to
#' the same units as \eqn{k_{elim}}.
#'
#' To convert to steady-state *blood* concentration, multiply by the
#' blood-to-plasma ratio.
#'
#' ### Oral route
#'
#' \deqn{C_{ss} = \frac{F_{gutabs} V_{dist}}{k_{elim} \tau}}
#'
#' ### Intravenous route
#'
#' \deqn{C_{ss} = \frac{1}{24 * \textrm{CL}_{tot}}}
#'
#' ## Half-life of elimination
#'
#' \deqn{\textrm{Halflife} = \frac{\log(2)}{k_{elim}}}
#'
#' ## Time of peak concentration
#'
#' For oral route:
#'
#' \deqn{\frac{\log \left( \frac{k_{gutabs}}{k_{elim}} \right)}{k_{gutabs} -
#' k_{elim}}}
#'
#' For intravenous route, time of peak concentration is always 0.
#'
#' ## Peak concentration
#'
#' Evaluate [cp_1comp()] at the time of peak concentration.
#'
#' ## AUC evaluated at infinite time
#'
#' Evaluate [auc_1comp()] at time = `Inf`.
#'
#' ## AUC evaluated at the time of the last observation
#'
#' Evaluate [auc_1comp()] at time = `tlast`.
#'
#'
#'
#'
#' @param pars A named vector of model parameters (e.g. from [coef.pk()]).
#' @param route Character: The route for which to compute TK stats. Currently
#'  only "oral" and "iv" are supported.
#' @param medium Character: the media (tissue) for which to compute TK stats.
#'  Currently only "blood" and "plasma" are supported.
#' @param dose Numeric: A dose for which to calculate TK stats.
#' @param time_unit Character: the units of time used for the parameters `par`.
#'  For example, if `par["kelim"]` is in units of 1/weeks, then `time_unit =
#'  "weeks"`. If `par["kelim"]` is in units of 1/hours, then `time_unit =
#'  "hours"`. This is used to calculate the steady-state plasma/blood
#'  concentration for long-term daily dosing of 1 mg/kg/day.
#' @param conc_unit Character: The units of concentration.
#' @param vol_unit Character: The units of dose.
#' @param ... Additional arguments not currently in use.
#' @return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css", "halflife", "tmax", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#' @export
#' @author John Wambaugh, Caroline Ring
tkstats_1comp <- function(pars,
                          route,
                          medium,
                          dose,
                          time_unit,
                          conc_unit,
                          vol_unit,
                          ...) {

  Fgutabs_Vdist <- Rblood2plasma <- NULL
  params <- fill_params_1comp(pars)

  # for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))


  CLtot <- kelim * Vdist

  CLtot_Fgutabs <- kelim / Fgutabs_Vdist

  # convert dose interval of (1/day) into time units
  # this is now standardized because time_units will always be hours
  dose_int <- 1 / 24

  Css <- dose * ifelse(route %in% "oral",
                      Fgutabs_Vdist / kelim / dose_int,
                      1 / (kelim * Vdist * dose_int)) *
    ifelse(medium %in% "blood",
           Rblood2plasma,
           1)

  halflife <- log(2) / kelim

  tmax <- ifelse(route %in% "oral",
                 log(kgutabs / kelim) / (kgutabs - kelim),
                 0)

  Cmax <- cp_1comp(params = pars,
  time = tmax,
  dose = dose,
  route = route,
  medium = medium,
  ...)

  AUC_inf <- auc_1comp(params = pars,
  time = Inf,
  dose = dose,
  route = route,
  medium = medium)

  return(data.frame(param_name = c("CLtot",
                                   "CLtot/Fgutabs",
                                   "Css",
                                   "halflife",
                                   "tmax",
                                   "Cmax",
                                   "AUC_infinity",
                                   "Vss",
                                   "Vss/Fgutabs"
                                   ),
                    param_value = c(CLtot,
                                    CLtot_Fgutabs,
                                    Css,
                                    halflife,
                                    tmax,
                                    Cmax,
                                    AUC_inf,
                                    Vdist,
                                    1 / (Fgutabs_Vdist)
                                    ),
                    param_units = c(paste0(vol_unit, "/", time_unit), # CLtot
                                    paste0(vol_unit, "/", time_unit), # CLtot/Fgutabs
                                    conc_unit, # Css
                                    time_unit, # halflife
                                    time_unit, # tmax
                                    conc_unit, # Cmax
                                    paste0(conc_unit, " * ", time_unit), # AUC_inf
                                    vol_unit,
                                    vol_unit)
                    ))

}
