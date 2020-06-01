logLikComp <- function(yout, ypred, M) {
    mu <- ypred %*% M
    return(sum(yout*log(mu)))
}


#' @title Permutation Test for Linear Independence Between Compositional Outcomes and Predictors
#'
#' @description Implements the loss function based permutation test as described in Fiksel et al. (2020)
#' for a test of linear indepence between compositional outcomes and predictors.
#'
#' @param yout A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param ypred A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param nperms The number of permutations. Default is 500.
#' @param init.seed The initial seed for the permutations. Default is 123.
#' @param accelerate A logical variable, indicating whether or not to use the Squarem algorithm for acceleration of the EM algorithm. Default is FALSE.
#'
#' @return The p-value for the indepence test
#' @export
codalm_indep_test <- function(yout, ypred, nperms = 500, init.seed = 123, accelerate = FALSE) {
    Nout <- nrow(yout)
    Npred <- nrow(ypred)
    if(Nout != Npred) {
        stop("Outcomes and predictors do not have the same number of observations")
    }
    for(i in 1:Nout) {
        yout[i,] <- yout[i,] / sum(yout[i,])
    }
    for(i in 1:Npred) {
        ypred[i,] <- ypred[i,] / sum(ypred[i,])
    }
    set.seed(init.seed)
    ### Get observed B
    B_obs <- codalm(yout, ypred, accelerate = accelerate)
    ### Get observed log-likelihood
    ll_full <- logLikComp(yout, ypred, B_obs)
    ### Get observed log-likelihood under null
    y_avg <- colMeans(yout)
    ll_null <- sum(sweep(yout, MARGIN = 2, log(y_avg), `*`))
    llr_obs <- ll_full - ll_null
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- sample(1:nrow(ypred))
        ypred_perm <- ypred[perm_index,]
        B_perm <- codalm(yout, ypred_perm, accelerate = accelerate)
        ll_perm <- logLikComp(yout, ypred_perm, B_perm)
        permut_stats[i] <- ll_perm - ll_null
    }
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}
