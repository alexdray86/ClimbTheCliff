#' MSE function used by CLIFF
#'
#' Callable with mse(tr, pr). Returns the mean square error between vector tr and pr.
#'
#' @param tr first input vector to be compared
#' @param pr second input vector to be compared
#' @export
mse <- function(tr, pr){ return( sum((tr - pr)^2) / length(tr) ) }
