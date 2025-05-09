#' Get scale_conc
#'
#' Extract concentration scale/transformation instructions from a [pk()] object
#'
#' @param obj A [pk()] object
#' @param ... Additional arguments not currently in use.
#' @return A `list`: `obj$scales$conc`
#' @author Caroline Ring
#' @export
get_scale_conc.pk <- function(obj, ...) {
  return(obj$scales$conc)
}
