#' Estimate parameter.
#' 
#' Estimate parameter from fitted model.
#'
#' @param object  fitted model
#' @param pars    names of parameters to estimate
#' @param ...     other parameters to pass
#' @export
estimate <- function(object, pars, ...) {
  UseMethod("estimate");
}

