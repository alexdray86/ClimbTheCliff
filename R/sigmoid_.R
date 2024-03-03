#' Sigmoid function with offset used by CLIFF
#' 
#' Apply the formula 1 / 1(1 + exp(-x + a)) 
#' Callable with sigmoid_(x, a)
#'
#' @param x input vector to be subjected to the sigmoid function
#' @param a shift to apply to the sigmoid function
#' @export
sigmoid_ <- function(x, a){
    return( 1 / (1 + exp(-x + a)))
}
