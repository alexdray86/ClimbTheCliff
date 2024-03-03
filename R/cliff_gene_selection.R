#' CLIFF gene selection
#'
#' Perform gene selection prior to CLIFF, using multiple rounds of lasso regression.
#'
#' @param bulk_mat matrix containing the input bulk matrix with gene expression
#' @param drug_d data.frame containing a column named 'auc' with AUC values of drug sensitivity
#' @param min.genes integer containing the minimum number of genes to consider
#' @export
cliff_gene_selection <- function(bulk_mat, drug_d, min.genes=10){
    stopifnot(dim(bulk_mat)[1] == length(drug_d))
    for(i in 1:100){
        la = 10^(-i)
        fit = glmnet(bulk_mat, y=drug_d, lambda=la)
        coefs_ = coef(fit)[-1]
        n_genes = sum(coefs_>0)
        if ( n_genes > min.genes ){
            return(colnames(bulk_mat)[coefs_ > 0])
        }
    }
}

