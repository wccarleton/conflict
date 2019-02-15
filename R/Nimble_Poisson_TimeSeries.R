poisTSCode <- nimbleCode({
   if(K > 1){
      A[1:K,1:K] <- inverse(diag(x=rep(1000,K)))
      multimu[1:K] <- rep(0,K)
      B[1:K] ~ dmnorm(mean=multimu[1:K], prec=A[1:K,1:K])
      mu[1] <- inprod(X[1,1:K], B[1:K])
   }else{
      B ~ dnorm(mean=0,sd=1000)
      mu[1] <- X[1]*B
   }
   sigma ~ dexp(0.1)
   lambda0 ~ dnorm(mean=0,sd=1000)
   rho ~ dnorm(0,sd=1000)
   lambda[1] ~ dnorm(rho * lambda0,sd=sigma)
   for (j in 2:T){
      if(K > 1){
         mu[j] <- inprod(X[j,1:K], B[1:K])
      }else{
         mu[j] <- X[j] * B
      }
      lambda[j] ~ dnorm(rho * lambda[j-1],sd=sigma)
   }
   for (j in 1:T){
      Y[j] ~ dpois(XBase[j]*exp(mu[j]+lambda[j]))
   }
})
#index <- 88:609 #379-900
#index <- index <- 309:609 #600-900
#index <- 1:706 #all
#index <- 59:609 #350–900
start <- 600
index <- c(1+(start-292):608)#292:900
#index <- c(1+(start-292):308)#start-600
#index <- 59:609 #350–900
#index <- 84:609 #375–900
#index <- 88:590 #379–881
#index <- 88:509 #379–800
#index <- 88:609 #379–900
#index <- 88:706 #379-end
#index <- c(72:597)#363-881
#index <- c(72:622)#363-913
Y <- MayaConflict_1y[index,2]
T <- length(Y)
XBase <- MayaConflict_1y[index,4]
SSTr = MayaConflict_1y[index,5]-mean(MayaConflict_1y[index,5])
SSTb = MayaConflict_1y[index,6]-mean(MayaConflict_1y[index,6])
d18O = MayaConflict_1y[index,7]-mean(MayaConflict_1y[index,7])
X <- cbind(SSTr,d18O,SSTr*d18O)
K <- ncol(X)

poisTSData <- list(Y=Y,
                  X=X,
                  XBase=XBase)

poisTSConsts <- list(T=T,
                     K=K)

poisTSInits <- list(lambda0=0,
                     sigma=1,
                     rho=0,
                     B=rep(0,K))

poisTSModel <- nimbleModel(code=poisTSCode,
                        data=poisTSData,
                        inits=poisTSInits,
                        constants=poisTSConsts)

#compile nimble model to C++ code—much faster runtime
C_poisTSModel <- compileNimble(poisTSModel, showCompilerOutput = FALSE)

#configure the MCMC
poisTSModel_conf <- configureMCMC(poisTSModel)

#select the variables that we want to monitor in the MCMC chain
poisTSModel_conf$addMonitors(c("lambda","logProb_Y"))

#build MCMC
poisTSModelMCMC <- buildMCMC(poisTSModel_conf,thin=1,enableWAIC=T)

#compile MCMC to C++—much faster
C_poisTSModelMCMC <- compileNimble(poisTSModelMCMC,project=poisTSModel)

#number of MCMC iterations
niter=200000

#set seed for replicability
set.seed(1000)

#call the C++ compiled MCMC model
C_poisTSModelMCMC$run(niter)

#save the MCMC chain (monitored variables) as a matrix
samples <- as.matrix(C_poisTSModelMCMC$mvSamples)
save(samples,file="../Results/MCMC_Chains/mcmc_model_3_interact.RData")
print(C_poisTSModelMCMC$calculateWAIC(nburnin=100000))
cols <- grep("log*",colnames(samples))
print(prod(colMeans(exp(samples[100000:200000,cols]))))
print(mean(samples[100000:200000,1]))
print(mean(samples[100000:200000,2]))
