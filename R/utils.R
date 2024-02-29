sigmoid_ <- function(x, a){
    return( 1 / (1 + exp(-x + a)))
}
mse <- function(tr, pr){ return( sum((tr - pr)^2) / length(tr) ) }
rmse <- function(tr, pr){ return( sqrt( sum((tr - pr)^2) / length(tr) ) ) }
