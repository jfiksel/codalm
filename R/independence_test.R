logLikComp <- function(y, x, M) {
    mu <- x %*% M
    return(sum(y*log(mu)))
}


#' @title Permutation Test for Linear Independence Between Compositional Outcomes and Predictors
#'
#' @description Implements the loss function based permutation test as described in Fiksel et al. (2020)
#' for a test of linear indepence between compositional outcomes and predictors.
#'
#' @param y A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param x A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param nperms The number of permutations. Default is 500.
#' @param init.seed The initial seed for the permutations. Default is 123.
#' @param accelerate A logical variable, indicating whether or not to use the
#' Squarem algorithm for acceleration of the EM algorithm. Default is TRUE.
#' @param parallel A logical variable, indicating whether or not to use a parallel
#' operation for computing the permutation statistics
#' @param ncpus Optional argument. When provided, is an integer giving the number
#' of clusters to be used in parallelization. Defaults to the number of cores, minus 1.
#'
#' @return The p-value for the indepence test
#' @export
#'
#' @import future
#' @import future.apply
#' @examples
#' require(gtools)
#' x <- rdirichlet(100, c(1, 1, 1))
#' y <- rdirichlet(100, c(1, 1, 1))
#' codalm_indep_test(y, x)
#'
#' require(ggtern)
#' data("WhiteCells", package = 'ggtern')
#' image <- subset(WhiteCells, Experiment == "ImageAnalysis")
#' image_mat <- as.matrix(image[,c("G", "L", "M")])
#' microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
#' microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])
#' x  <- image_mat  / rowSums(image_mat)
#' y <- microscopic_mat / rowSums(microscopic_mat)
codalm_indep_test <- function(y, x, nperms = 500, init.seed = 123, accelerate = TRUE,
                              parallel = FALSE, ncpus = NULL) {
    Nout <- nrow(y)
    Npred <- nrow(x)
    if(Nout != Npred) {
        stop("Outcomes and predictors do not have the same number of observations")
    }
    for(i in 1:Nout) {
        y[i,] <- y[i,] / sum(y[i,])
    }
    for(i in 1:Npred) {
        x[i,] <- x[i,] / sum(x[i,])
    }
    ### Get observed B
    B_obs <- codalm(y, x, accelerate = accelerate)
    ### Get observed log-likelihood
    ll_full <- logLikComp(y, x, B_obs)
    ### Get observed log-likelihood under null
    y_avg <- colMeans(y)
    ll_null <- sum(sweep(y, MARGIN = 2, log(y_avg), `*`))
    llr_obs <- ll_full - ll_null
    ### Get do permutation test
    RNGkind("L'Ecuyer-CMRG")
    if(parallel == FALSE) {
        plan(sequential)
    } else {
        if(is.null(ncpus)) {
            nworkers <- availableCores() - 1
        } else {
            nworkers <- min(availableCores() - 1, ncpus)
        }
        plan(multiprocess, workers = nworkers)

    }
    permut_stats <- future_sapply(1:nperms, function(i) {
        perm_index <- sample(1:nrow(x))
        x_perm <- x[perm_index,]
        B_perm <- codalm(y, x_perm, accelerate = accelerate)
        ll_perm <- logLikComp(y, x_perm, B_perm)
        return(ll_perm - ll_null)
    }, future.seed = init.seed)
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}
