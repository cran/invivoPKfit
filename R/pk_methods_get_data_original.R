#' Get data_original
#'
#' @param obj A [pk()] object
#' @param ... Additional arguments. Not currently in use.
#' @return A `data.frame` -- the `data_original` element of `obj`
#' @export
#' @author Caroline Ring
get_data_original.pk <- function(obj, ...) {
obj$data_original
}
