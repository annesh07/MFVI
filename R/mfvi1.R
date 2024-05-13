#' Mean-Field Variational Inference-I
#'
#' This function performs the mean-field variational inference for a DP mixture model, where the alpha parameter from the DP prior is considered varied (a Gamma prior) and the variance of the observations is considered fixed
#'
#' @param D Number of Dimension
#' @param N Number of Observations
#' @param C0 Covariance matrix of the Observed data
#' @param T0 Number of clusters in the variational distribution
#' @param s1 Shape parameter of the prior Gamma distribution for alpha
#' @param s2 Rate parameter of the prior Gamma distribution for alpha
#' @param X Observations, N x D matrix
#' @param L20 Precision parameter of the prior multivariate normal distribution for eta
#' @param Mu0 Mean parameter of the prior multivariate normal distribution for eta
#' @param W1 Shape parameter of the posterior gamma distribution for alpha
#' @param W2 Rate parameter of the posterior gamma distribution for alpha
#' @param V1 1st shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param V2 2nd shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param L1 Mean parameters of the posterior multivariate normal distribution for eta's, T0 x D matrix
#' @param L2 Precision parameters of the posterior multivariate normal distribution for eta's, T0 x 1 matrix
#' @param Plog Posterior probabilities of the latent allocations
#' @param maxit Maximum number of iterations for the variational inference
#'
#' @return the posterior expected value of alpha and the number of clusters on the basis of latent allocation probabilities
#' @export
#'
#' @examples mfvi1(D=2,N=100,C0=diag(2),T0=20,s1=0.001,s2=0.001,X=matrix(1,nrow=100, ncol=2),L20=0.001,Mu0=matrix(c(rep(0,2)), nrow=1),W1=0.001,W2=0.001,V1=matrix(0.0001, nrow = 20, ncol = 1),V2=matrix(0.0001, nrow = 20, ncol = 1),L1=matrix(0.0001, nrow = 20, ncol = 2),L2=matrix(0.0001, nrow = 20, ncol = 1),Plog=matrix(-3,nrow=100,ncol=20),maxit=1000)
mfvi1 <- function(D,N,C0,T0,s1,s2,X,L20,Mu0,W1,W2,V1,V2,L1,L2,Plog,maxit){
  C00 <- diag(D)/L20
  Mu00 <- Mu0%*%solve(C00)
  source("C:/Users/ap15/Documents/MFVI/R/elbo1.R")

  f <- list()
  f[[1]] = ELBO1(D,N,C0,T0,s1,s2,X,L20,Mu0,W1,W2,V1,V2,L1,L2,Plog)
  for (m in 1:maxit){
    W1 = s1 + T0
    W2 = s2 - sum(digamma(V2)-digamma(V1+V2))

    for (j in 1:T0){
      V1[j,1] = 1 + sum(exp(Plog[,j]))
      j = j+1
    }

    for (j0 in 1:T0-1){
      V2[j0,1] = W1/W2 + sum(exp(Plog[,(j0+1):T0]))
      j0 = j0+1
    }
    V2[T0,1] = W1/W2

    for (j in 1:T0){
      L1[j,] = Mu00 + t(exp(Plog[,j, drop=FALSE]))%*%X
      j = j+1
    }
    for (j in 1:T0){
      L2[j,1] = L20 + sum(exp(Plog[,j]))
      j = j+1
    }
    for (i in 1:N){
      Plog[i,1] = digamma(V1[1,1]) - digamma(V1[1,1]+V2[1,1]) + (L1[1,, drop=FALSE]/L2[1,1])%*%t(X[i,, drop=FALSE]) - 0.5*((L1[1,, drop=FALSE]/L2[1,1])%*%(t(L1[1,, drop=FALSE])/L2[1,1])+D/L2[1,1])
      for (j00 in 2:T0){
        Plog[i,j00] = digamma(V1[j00,1]) - digamma(V1[j00,1]+V2[j00,1]) + sum(digamma(V2[1:(j00-1),1]) - digamma(V1[1:(j00-1),1]+V2[1:(j00-1),1])) + (L1[j00,, drop=FALSE]/L2[j00,1])%*%t(X[i,, drop=FALSE]) - 0.5*((L1[j00,, drop=FALSE]/L2[j00,1])%*%(t(L1[j00,, drop=FALSE])/L2[j00,1])+D/L2[j00,1])
        j00 = j00+1
      }
      p0 = max(Plog[i,])
      Plog[i,] = Plog[i,]-p0-log(sum(exp(Plog[i,]-p0)))
      i = i+1
    }

    f[[m+1]] = ELBO1(D,N,C0,T0,s1,s2,X,L20,Mu0,W1,W2,V1,V2,L1,L2,Plog)
    if (abs(sum(f[[m]])-sum(f[[m+1]]))<0.000001){
      break
    }
    cat("outer loop: ", m,"\n", sep="")
    print(f[[m+1]])
    cat('\n')
  }
  alpha <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clustnum <- length(unique(clustering))
  list0 <- list("alpha"=alpha, "Clusters"=clustnum)
  return(list0)
}
