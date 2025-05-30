#' Helper function to convert time units
#'
#' Convert a vector of times between units
#'
#'
#' @param x Numeric: one or more time values to be converted.
#' @param from Character vector: `x` is currently in these units. Must be units
#'   understood by `lubridate::duration()`, i.e. `"seconds"`, `"hours"`,
#'   `"days"`, `"weeks"`, `"months"`, `"years"`, `"milliseconds"`,
#'   `"microseconds"`, `"nanoseconds"`, and/or `"picoseconds"`. Default value is
#'   `"hours"`.
#' @param to Character vector: `x` will be converted to these units. Must be
#'   either `"auto"`, `"identity"`, or units understood by
#'   `lubridate::duration()`, i.e. `"seconds"`, `"hours"`, `"days"`, `"weeks"`,
#'   `"months"`, `"years"`, `"milliseconds"`, `"microseconds"`, `"nanoseconds"`,
#'   and/or `"picoseconds"`. Default value is `"identity"`.  If `"identity"`,
#'   then `x` will be returned unchanged. If `"auto"`, then units will be
#'   automatically chosen that make the midpoint of `x` (or its inverse, if
#'   `inverse = TRUE`) as close to an order of magnitude of 10 as possible (see
#'   [auto_units()]).
#' @param inverse Logical: TRUE if `x` is in units of *inverse* time (e.g.
#'   1/hour, 1/day); FALSE if `x` is in units of time (e.g. hours, days).
#'   Default value is FALSE.
#' @return A numeric vector the same length as `x`, converted from the units in
#'   `from` to the units in `to`.
#'
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado

convert_time <- function(x,
                         from = "hours",
                         to = "identity",
                         inverse = FALSE) {

  period_units <- time_units
  time_conversions <- time_conversions

  if (to %in% "auto") {
    # select units automatically
    to <- auto_units(y = x, from = from)
  }


  if (to != "identity") {
    if (isTRUE(inverse)) {
      y <- 1 / x
    } else {
      y <- x
    }

    if (to %in% period_units) {
      # test with convert_time(c(60, 15, 30, 120), from = "minutes", to = "hours")
      conversion_table <- time_conversions[
        to == time_conversions$TimeTo &
          unique(from) == time_conversions$TimeFrom, "conversion"
      ]

      y_new <- y * conversion_table

    } else {
      y_new <- NA_real_
      message(
          "invivopkfit::convert_time(): Allowable values for 'to' are",
          toString(period_units),
          ", but 'to' =",
          to,
          ". Returning NA_real_ for converted times."
      )
    }

    if (isTRUE(inverse)) {
      yout <- 1 / y_new
    } else {
      yout <- y_new
    }
  } else {
    # if to == "identity", then just return the input as-is
    yout <- x
  }

  return(yout)
}
