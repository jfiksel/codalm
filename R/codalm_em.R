softmax <- function(x){exp(x)/sum(exp(x))}


compreg.em <- function(pars, A, G, D1, D2) {
    ### Note that pars should be a C^2 vector representing M
    M <- matrix(pars, ncol = D1, nrow = D2, byrow = FALSE)
    Mnew <- matrix(NA, ncol = D1, nrow = D2)
    ### First compute generic weights for E[Z_i | d = j]
    weights_array <- array(NA, dim = c(nrow(A), D2, D1))
    ### W_{rij} = weights_array[r,i,j]
    for(j in 1:D1) {
        weights <- sweep(G, MARGIN=2, M[,j], `*`)
        weights <- weights/rowSums(weights)
        weights_array[,,j] <- weights
    }
    for(i in 1:D2) {
        for(j in 1:D1) {
            w_ij <- weights_array[,i,j]
            Mnew[i,j] <- sum(A[,j] * w_ij)
        }
    }
    Mnew[Mnew<=0] <- 1e-8
    Mnew <- Mnew/rowSums(Mnew)
    parsnew <- as.vector(Mnew)
    return(parsnew)
}


compreg.loglik <- function(pars, A, G, D1, D2) {
    N <- nrow(A)
    M <- matrix(pars, ncol = D1, nrow = D2, byrow = FALSE)
    mu <- G %*% M
    loglik <- 0
    for(r in 1:N) {
        for(j in 1:D1) {
            loglik <- loglik + A[r,j] * log(mu[r,j] + .001)
        }
    }
    return(-loglik)
}


#' @title Transformation-free Linear Regression for Compositional Outcomes and Predictors
#'
#' @description Implements the expectation-maximization (EM) algorithm as described
#' in Fiksel et al. (2020) for transformation-free linear regression for
#' compositional outcomes and predictors.
#'
#' @param y A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param x A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param accelerate A logical variable, indicating whether or not to use the
#' Squarem algorithm for acceleration of the EM algorithm. Default is TRUE.
#'
#' @references \url{https://arxiv.org/abs/2004.07881}
#'
#' @return A \eqn{D_s} x \eqn{D_r} compositional coefficient matrix, where
#' \eqn{D_s} and \eqn{D_r}  are the dimensions of the compositional predictor
#' and outcome, respectively
#' @export
#' @importFrom SQUAREM squarem fpiter
#' @examples
#' require(ggtern)
#' data("WhiteCells", package = 'ggtern')
#' image <- subset(WhiteCells, Experiment == "ImageAnalysis")
#' image_mat <- as.matrix(image[,c("G", "L", "M")])
#' microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
#' microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])
#' x <- image_mat  / rowSums(image_mat)
#' y <- microscopic_mat / rowSums(microscopic_mat)
#' codalm(y, x)
codalm <- function(y, x, accelerate = TRUE) {
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
    D1 <- ncol(y)
    D2 <- ncol(x)
    par0 <- rep(1/D1, D1*D2)
    if(accelerate) {
        em_output <- squarem(par = par0, fixptfn = compreg.em,
                             objfn = compreg.loglik, control = list(tol = 1.e-8),
                             A = y, G = x, D1 = D1, D2 = D2)
    } else {
        em_output <- fpiter(par = par0, fixptfn = compreg.em,
                            objfn = compreg.loglik, control = list(tol = 1.e-8),
                            A = y, G = x, D1 = D1, D2 = D2)
    }
    B_est <- em_output$par
    dim(B_est) <- c(D2, D1)
    return(B_est)
}

codalm_boot_fn <- function(y, x, indices, accelerate) {
    yboot <- y[indices,]
    xboot <- x[indices,]
    B_est <- codalm(yboot, xboot, accelerate = accelerate)
    return(B_est)
}

#' @title Bootstrap Confidence Intervals Linear Regression for Compositional Outcomes and Predictors
#'
#' @description Implements percentile based bootstrapping to estimate the confidence intervals
#' for the regression coefficients when doing linear regression for compositional outcomes
#' and predictors
#'
#' @param y A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param x A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param accelerate A logical variable, indicating whether or not to use the
#' Squarem algorithm for acceleration of the EM algorithm. Default is TRUE
#' @param nboot The number of bootstrap repetitions to use. Default is 500
#' @param conf A scalar between 0 and 1 containing the confidence level of the required intervals.
#' Default is .95.
#' @param parallel A logical variable, indicating whether or not to use a parallel
#' operation for computing the permutation statistics
#' @param ncpus Optional argument. When provided, is an integer giving the number
#' of clusters to be used in parallelization. Defaults to the number of cores, minus 1.
#' @param strategy Optional argument. When provided, this will be the evaluation function
#' (or name of it) to use for parallel computation (if parallel = TRUE). Otherwise,
#' if parallel = TRUE, then this will default to multisession. See \code{\link[future]{plan}}.
#' @param init.seed The initial seed for the permutations. Default is 123.
#'
#' @return A list, with \code{ci_L} and \code{ci_U}, giving the lower and upper bounds
#' of each element of the B matrix
#' @export
#'
#' @import future
#' @import future.apply
#' @importFrom stats quantile
#' @examples
#' require(ggtern)
#' data("WhiteCells", package = 'ggtern')
#' image <- subset(WhiteCells, Experiment == "ImageAnalysis")
#' image_mat <- as.matrix(image[,c("G", "L", "M")])
#' microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
#' microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])
#' x <- image_mat  / rowSums(image_mat)
#' y <- microscopic_mat / rowSums(microscopic_mat)
#' codalm_ci(y, x, nboot = 50, conf = .95)
codalm_ci <- function(y, x, accelerate = TRUE, nboot = 500, conf = .95,
                      parallel = FALSE, ncpus = NULL,
                      strategy = NULL, init.seed = 123) {
    ### Get future strategy
    if(parallel == FALSE) {
        strategy <- 'sequential'
    } else  {
        strategy <- ifelse(is.null(strategy), 'multisession', strategy)
        ### get number of workers
        if(is.null(ncpus)) {
            nworkers <- availableCores() - 1
        } else {
            nworkers <- min(availableCores() - 1, ncpus)
        }
        ### make sure nworkers isnt 0
        nworkers <- ifelse(nworkers <= 0, 1, nworkers)
    }
    ### set plan and make sure future resets at end of function
    if(strategy == 'sequential') {
        oplan <- plan(strategy)
    } else {
        oplan <- plan(strategy, workers = nworkers)
    }
    ### make sure future resets at end of function
    on.exit(plan(oplan), add = TRUE)
    Nout <- nrow(y)
    Npred <- nrow(x)
    if(Nout != Npred) {
        stop("Outcomes and predictors do not have the same number of observations")
    }
    ymat <- cbind(y, x)
    D1 <- ncol(y)
    D2 <- ncol(x)
    ### Do permutation test
    RNGkind("L'Ecuyer-CMRG")
    bootstraps <- future_lapply(1:nboot, function(i) {
        boot_index <- sample(1:nrow(y), replace = TRUE)
        yboot <- y[boot_index,]
        xboot <- x[boot_index,]
        B_est <- codalm(yboot, xboot, accelerate = accelerate)
        return(B_est)
    }, future.seed = init.seed)
    prob_L <- (1 - conf)/2
    prob_U <- conf + prob_L
    quantiles <- apply(simplify2array(bootstraps), 1:2, quantile, prob = c(prob_L, prob_U))
    ci_L <- quantiles[1,,]
    ci_U <- quantiles[2,,]
    return(list(ci_L = ci_L, ci_U = ci_U))
}
