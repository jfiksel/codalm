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
#' @param yout A matrix of compositional outcomes. Each row is an observation, and must sum to 1
#' @param ypred A matrix of compositional predictors. Each row is an observation, and must sum to 1
#' @param accelerate A logical variable, indicating whether or not to use the Squarem algorithm for acceleration of the EM algorithm. Default is FALSE.
#'
#' @references \url{https://arxiv.org/abs/2004.07881}
#'
#' @return A \eqn{D_s} x \eqn{D_r} compositional coefficient matrix, where \eqn{D_s} and \eqn{D_r}  are the dimensions of the compositional predictor and outcome, respectively
#' @export
#' @importFrom SQUAREM squarem fpiter
codalm <- function(yout, ypred, accelerate = FALSE) {
    D1 <- ncol(yout)
    D2 <- ncol(ypred)
    par0 <- rep(1/D1, D1*D2)
    if(accelerate) {
        em_output <- squarem(par = par0, fixptfn = compreg.em,
                             objfn = compreg.loglik, control = list(tol = 1.e-8),
                             A = yout, G = ypred, D1 = D1, D2 = D2)
    } else {
        em_output <- fpiter(par = par0, fixptfn = compreg.em,
                            objfn = compreg.loglik, control = list(tol = 1.e-8),
                            A = yout, G = ypred, D1 = D1, D2 = D2)
    }
    B_est <- em_output$par
    dim(B_est) <- c(D2, D1)
    return(B_est)
}
