#' RMSE function used by CLIFF
#'
#' Callable with rmse(tr, pr). Returns the root mean square error between vector tr and pr.
#'
#' @param tr first input vector to be compared
#' @param pr second input vector to be compared
#' @export
rmse <- function(tr, pr){ return( sqrt( sum((tr - pr)^2) / length(tr) ) ) }
