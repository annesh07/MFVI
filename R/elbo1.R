#' ELBO function for DPMM with varied alpha and fixed sigma
#'
#' Computes the ELBO for a DPMM with varied alpha and fixed sigma
#'
#' @param d Number of Dimensions
#' @param n Number of Observations
#' @param C0 Covariance matrix of the Observed data
#' @param T0 Number of clusters in the variational distribution
#' @param shape_alpha Shape parameter of the prior Gamma distribution for alpha
#' @param rate_alpha Rate parameter of the prior Gamma distribution for alpha
#' @param X Observation matrix, n x d
#' @param L20 Precision parameter of the prior multivariate normal distribution for eta
#' @param Mu0 Mean parameter of the prior multivariate normal distribution for eta
#' @param W1 Shape parameter of the posterior gamma distribution for alpha
#' @param W2 Rate parameter of the posterior gamma distribution for alpha
#' @param V1 1st shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param V2 2nd shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param L1 Mean parameters of the posterior multivariate normal distribution for eta's, T0 x d matrix
#' @param L2 Precision parameters of the posterior multivariate normal distribution for eta's, T0 x 1 matrix
#' @param Plog Posterior probabilities of the latent allocations
#'
#' @return the value of the ELBO
#' @export
#'
#' @examples
#' ELBO1(d=2, n=100,C0=diag(2),T0=20,shape_alpha=0.001,rate_alpha=0.001,X=matrix(1,nrow=100, ncol=2),L20=0.001,Mu0=matrix(c(rep(0,2)), nrow=1),W1=0.001,W2=0.001,V1=matrix(0.0001, nrow = 20, ncol = 1),V2=matrix(0.0001, nrow = 20, ncol = 1),L1=matrix(0.0001, nrow = 20, ncol = 2),L2=matrix(0.0001, nrow = 20, ncol = 1),Plog=matrix(-3,nrow=100,ncol=20))
ELBO1 <- function(d, n, C0,T0,shape_alpha,rate_alpha,X,L20,Mu0,W1,W2,V1,V2,L1,L2,Plog){
  C00 <- diag(d)/L20
  V1 = matrix(V1, nrow = T0, ncol = 1)
  V2 = matrix(V2, nrow = T0, ncol = 1)
  L1 = matrix(L1, nrow = T0, ncol = d)
  L2 = matrix(L2, nrow = T0, ncol = 1)
  Plog = matrix(Plog, nrow = n, ncol = T0)

  #the alpha
  e0 <- shape_alpha*log(rate_alpha) - lgamma(shape_alpha) + (shape_alpha-1)*(-log(W2) + digamma(W1)) - rate_alpha*W1/W2
  #the Vi's
  e1 <- -T0*lbeta(1, W1/W2)+(W1/W2-1)*sum(digamma(V2)-digamma(V1+V2))
  #the etai's
  e20 = c()
  inv_C00 = solve(C00)
  for (i in 1:T0){
    e20[i] = -d/2*log(2*pi) + d*0.5*log(L20) - 0.5*(((L1[i,, drop=FALSE]/L2[i,1])%*%inv_C00%*%(t(L1[i,, drop=FALSE])/L2[i,1])) + sum(diag(inv_C00/L2[i,1]))) + Mu0%*%inv_C00%*%(t(L1[i,, drop=FALSE])/L2[i,1]) - 0.5*Mu0%*%inv_C00%*%t(Mu0)
  }
  e2 = sum(e20)
  #the zn's
  P0 = matrix(0, nrow = n, ncol = T0)
  for (i in 1:n){
    for (j in 1:(T0-1)){
      P0[i,j] = sum(exp(Plog[i,(j+1):T0]))
    }
    P0[i,T0] = 0
  }
  e3 = sum(sweep(P0,2,digamma(V2)-digamma(V2+V1),"*")+sweep(exp(Plog),2,digamma(V1)-digamma(V2+V1),"*"))
  #the xn's
  e40 = matrix(NA, nrow = n, ncol = T0)
  inv_C0 = solve(C0)
  for (i in 1:n){
    for (j in 1:T0){
      e40[i,j] = exp(Plog[i,j])*(X[i,, drop=FALSE]%*%inv_C0%*%t(L1[j,, drop=FALSE]/L2[j,1])-0.5*X[i,, drop=FALSE]%*%inv_C0%*%t(X[i,, drop=FALSE])-0.5*(d*log(2*pi) + determinant(C0, logarithm=TRUE)$modulus + (L1[j,, drop=FALSE]/L2[j,1])%*%inv_C0%*%t(L1[j,, drop=FALSE]/L2[j,1])) + sum(diag(inv_C0/L2[j,1])))
    }
  }
  e4 = sum(e40)
  #the variational parameters
  e5 = sum(-lbeta(V1,V2) + (V1-1)*(digamma(V1)-digamma(V2+V1)) + (V2-1)*(digamma(V2)-digamma(V2+V1)))
  e50 = sum(exp(Plog)*Plog)
  e51 = c()
  for (i in 1:T0){
    e51[i] = -d/2*log(2*pi) + d*0.5*log(L2[i,1]) - 0.5*(((L1[i,, drop=FALSE]/L2[i,1])%*%(t(L1[i,, drop=FALSE]/L2[i,1])))*L2[i,1]+d) + (0.5*L2[i,1])*((L1[i,, drop=FALSE]/L2[i,1])%*%t(L1[i,, drop=FALSE]/L2[i,1]))
  }
  e52 = W1*log(W2) - lgamma(W1) +(W1-1)*(-log(W2) + digamma(W1)) - W1
  e5 = e5 + e50 + sum(e51) + e52
  #browser()
  #e = e1 + e2 + e3 + e4 - e5
  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
}
