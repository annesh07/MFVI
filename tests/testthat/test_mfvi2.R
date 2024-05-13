source("C:/Users/ap15/Documents/MFVI/R/mfvi2.R")
test_that("Number of clusters from latent allocation:",{
  expect_equal(mfvi2(d=2,n=100,T0=20,shape_alpha=0.001,rate_alpha=0.001,s01=0.001,s02=0.001,X=matrix(1,nrow=100, ncol=2),L20=0.001,Mu0=matrix(c(rep(0,2)), nrow=1),W1=0.001,W2=0.001,G1=0.001,G2=0.001,V1=matrix(0.0001, nrow = 20, ncol = 1),V2=matrix(0.0001, nrow = 20, ncol = 1),L1=matrix(0.0001, nrow = 20, ncol = 2),L2=matrix(0.0001, nrow = 20, ncol = 1),Plog=matrix(-3,nrow=100,ncol=20),maxit=1000)$Clusters,1)
  })
