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
#' @description Implements the EM algorithm as described in Fiksel et al. (2020)
#' for transformation-free linear regression for compositional outcomes and predictors.
#'
#' @param y A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param x A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param accelerate A logical variable, indicating whether or not to use the
#' Squarem #algorithm for acceleration of the EM algorithm. Default is TRUE.
#'
#' @references \url{https://arxiv.org/abs/2004.07881}
#'
#' @return A \eqn{D_s} x \eqn{D_r} compositional coefficient matrix, where \eqn{D_s} and \eqn{D_r}  are the dimensions of the compositional predictor and outcome, respectively
#' @export
#' @importFrom SQUAREM squarem fpiter
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

codalm_boot_fn <- function(ymat, indices, D1, D2, accelerate) {
    y <- ymat[indices, 1:D1]
    x <- ymat[indices, (D1+1):(D1+D2)]
    B_est <- codalm(y, x, accelerate = accelerate)
    return(as.vector(B_est))
}

#' @title Bootstrap Confidence Intervals Linear Regression for Compositional Outcomes and Predictors
#'
#' @description Implements the EM algorithm as described in Fiksel et al. (2020)
#' for transformation-free linear regression for compositional outcomes and predictors.
#'
#' @param y A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param x A matrix of compositional predictors. Each row is an observation, and must sum to 1.
#' If any rows do not sum to 1, they will be renormalized
#' @param accelerate A logical variable, indicating whether or not to use the
#' Squarem algorithm for acceleration of the EM algorithm. Default is TRUE
#' @param R The number of bootstrap reptitions to use. Default is 500
#' @param ci_type A character string with the type of bootstrap interval to be calculated. One of
#' "norm", "perc", "basic", or "bca". See the documentation for
#' \code{\link[boot]{boot.ci}} for more information. Default is "perc".
#' @param conf A scalar between 0 and 1 containing the confidence level of the required intervals.
#' Default is .95.
#' @param ... Additional arguments to pass to the \code{\link[boot]{boot}} function.
#' If you want to use parallelize the bootstrapping, you can pass in the
#' appropriate arguments here.
#'
#' @return A list, with \code{ci_L} and
#' @export
#'
#' @importFrom boot boot boot.ci
coda_lm_ci <- function(y, x, accelerate = TRUE, R = 500, ci_type = "perc", conf = .95, ...) {
    Nout <- nrow(y)
    Npred <- nrow(x)
    if(Nout != Npred) {
        stop("Outcomes and predictors do not have the same number of observations")
    }
    ymat <- cbind(y, x)
    D1 <- ncol(y)
    D2 <- ncol(x)
    bootstraps <- boot(ymat, codalm_boot_fn, R = R, D1 = D1, D2 = D2, accelerate = accelerate, ...)
    ci_mat <- do.call(rbind, lapply(1:(D1*D2), function(i) {
        ci <- boot.ci(bootstraps, index = i, conf = conf, type = c("norm", "basic", "perc", "bca"))
        if(ci_type == "norm") {
            return(ci$normal[,2:3])
        } else if (ci_type == "basic") {
            return(ci$basic[,4:5])
        } else if (ci_type == "perc") {
            return(ci$percent[,4:5])
        } else {
            return(ci$bca[,4:5])
        }
    }))
    ci_L <- as.vector(ci_mat[,1])
    dim(ci_L) <- c(D2, D1)
    ci_U<- as.vector(ci_mat[,2])
    dim(ci_U) <- c(D2, D1)
    return(list(ci_L = ci_L, ci_U = ci_U))
}
