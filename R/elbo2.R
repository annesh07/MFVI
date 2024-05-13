#' ELBO function for DPMM with varied alpha and varied sigma
#'
#' This function computes the ELBO for a DPMM with varied alpha and varied sigma
#'
#' @param D Number of Dimensions
#' @param N Number of Observations
#' @param T0 Number of clusters in the variational distribution
#' @param s1 Shape parameter of the prior Gamma distribution for alpha
#' @param s2 Rate parameter of the prior Gamma distribution for alpha
#' @param s01 Shape parameter of the prior gamma distribution for sigma
#' @param s02 Rate parameter of the prior gamma distribution for sigma
#' @param X Observations, N x D matrix
#' @param L20 Precision parameter of the prior multivariate normal distribution for eta
#' @param Mu0 Mean parameter of the prior multivariate normal distribution for eta
#' @param W1 Shape parameter of the posterior gamma distribution for alpha
#' @param W2 Rate parameter of the posterior gamma distribution for alpha
#' @param G1 Shape parameter of the posterior gamma distribution for sigma
#' @param G2 Rate parameter of the posterior gamma distribution for sigma
#' @param V1 1st shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param V2 2nd shape parameters of the posterior beta distribution for v's, T0 x 1 matrix
#' @param L1 Mean parameters of the posterior multivariate normal distribution for eta's, T0 x D matrix
#' @param L2 Precision parameters of the posterior multivariate normal distribution for eta's, T0 x 1 matrix
#' @param Plog Posterior probabilities of the latent allocations
#'
#' @return the value of the ELBO
#' @export
#'
#' @examples
#' ELBO2(D=2,N=100,T0=20,s1=0.001,s2=0.001,s01=0.001,s02=0.001,X=matrix(1,nrow=100, ncol=2),L20=0.001,Mu0=matrix(c(rep(0,2)), nrow=1),W1=0.001,W2=0.001,G1=0.001,G2=0.002,V1=matrix(0.0001, nrow = 20, ncol = 1),V2=matrix(0.0001, nrow = 20, ncol = 1),L1=matrix(0.0001, nrow = 20, ncol = 2),L2=matrix(0.0001, nrow = 20, ncol = 1),Plog=matrix(-3,nrow=100,ncol=20))
ELBO2 <- function(D,N,T0,s1,s2,s01,s02,X,L20,Mu0,W1,W2,G1,G2,V1,V2,L1,L2,Plog){
  C00 <- diag(D)/L20
  W1 = W1
  W2 = W2
  G1 = G1
  G2 = G2
  V1 = matrix(V1, nrow = T0, ncol = 1)
  V2 = matrix(V2, nrow = T0, ncol = 1)
  L1 = matrix(L1, nrow = T0, ncol = D)
  L2 = matrix(L2, nrow = T0, ncol = 1)
  Plog = matrix(Plog, nrow = N, ncol = T0)

  #the alpha
  e0 = s1*log(s2) - lgamma(s1) + (s1-1)*(-log(W2) + digamma(W1)) - s2*W1/W2
  #the Vi's
  e1 = -T0*lbeta(1, W1/W2)+(W1/W2-1)*sum(digamma(V2)-digamma(V1+V2))
  #the etai's
  e20 = c()
  inv_C00 = solve(C00)
  for (i in 1:T0){
    e20[i] = -D/2*log(2*pi) + D*0.5*log(L20) - 0.5*(((L1[i,, drop=FALSE]/L2[i,1])%*%inv_C00%*%(t(L1[i,, drop=FALSE])/L2[i,1])) + sum(diag(inv_C00/L2[i,1]))) + Mu0%*%inv_C00%*%(t(L1[i,, drop=FALSE])/L2[i,1]) - 0.5*Mu0%*%inv_C00%*%t(Mu0)
  }
  e2 = sum(e20)
  #the zn's
  P0 = matrix(0, nrow = N, ncol = T0)
  for (i in 1:N){
    for (j in 1:(T0-1)){
      P0[i,j] = sum(exp(Plog[i,(j+1):T0]))
    }
    P0[i,T0] = 0
  }
  e3 = sum(sweep(P0,2,digamma(V2)-digamma(V2+V1),"*")+sweep(exp(Plog),2,digamma(V1)-digamma(V2+V1),"*"))
  #the covariance
  e4 = s01*log(s02) - lgamma(s01) + (s01-1)*(-log(G2) + digamma(G1)) - s02*G1/G2
  #the xn's
  e50 = matrix(NA, nrow = N, ncol = T0)
  #inv_C0 = solve(C0)
  for (i in 1:N){
    for (j in 1:T0){
      e50[i,j] = exp(Plog[i,j])*((G1/G2)*X[i,, drop=FALSE]%*%t(L1[j,, drop=FALSE]/L2[j,1]) - 0.5*(G1/G2)*X[i,, drop=FALSE]%*%t(X[i,, drop=FALSE]) - 0.5*D*log(2*pi) + 0.5*D*log(G1/G2) - 0.5*(G1/G2)*((L1[j,, drop=FALSE]/L2[j,1])%*%t(L1[j,, drop=FALSE]/L2[j,1]) + D/L2[j,1]))
    }
  }
  e5 = sum(e50)
  #the priors
  e6 = sum(-lbeta(V1,V2) + (V1-1)*(digamma(V1)-digamma(V2+V1)) + (V2-1)*(digamma(V2)-digamma(V2+V1)))
  e60 = sum(exp(Plog)*Plog)
  e61 = c()
  for (i in 1:T0){
    e61[i] = -D/2*log(2*pi) + D*0.5*log(L2[i,1]) - 0.5*(((L1[i,, drop=FALSE]/L2[i,1])%*%(t(L1[i,, drop=FALSE]/L2[i,1])))*L2[i,1]+D) + (0.5*L2[i,1])*((L1[i,, drop=FALSE]/L2[i,1])%*%t(L1[i,, drop=FALSE]/L2[i,1]))
  }
  e62 = W1*log(W2) - lgamma(W1) +(W1-1)*(-log(W2) + digamma(W1)) - W1
  e63 = G1*log(G2) - lgamma(G1) +(G1-1)*(-log(G2) + digamma(G1)) - G1
  e6 = e6 + e60 + sum(e61) + e62 + e63
  #browser()
  #e = e1 + e2 + e3 + e4 - e5
  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "e5"=e5, "me6"=-e6))
}
