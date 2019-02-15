poisTSCode <- nimbleCode({
   sigma ~ dexp(0.1)
   lambda0 ~ dnorm(mean=0,sd=1000)
   rho ~ dnorm(0,sd=1000)
   lambda[1] ~ dnorm(rho * lambda0,sd=sigma)
   for (j in 2:T){
      lambda[j] ~ dnorm(rho * lambda[j-1],sd=sigma)
   }
   for (j in 1:T){
      Y[j] ~ dpois(XBase[j]*exp(lambda[j]))
   }
})
#index <- 88:609 #379-900
#index <- 309:609 #600-900
#index <- 1:706 #all
#index <- c(88:590)
#index <- c(1:609)#292-900
start <- 292
#index <- c(1+(start-292):87)#start-379
#index <- c(1+(start-292):308)#start-600
index <- c(1+(start-292):608)#start-900
#index <- 88:609 #379-900
#index <- 88:509 #379–800
#index <- 88:706 #379-end
#index <- 59:609 #350–900
#index <- 88:590 #379–881
#index <- c(72:597)#363-881
#index <- c(72:622)#363-913
Y <- MayaConflict_1y[index,2]
T <- length(Y)
XBase <- MayaConflict_1y[index,4]

poisTSData <- list(Y=Y,
                  XBase=XBase)

poisTSConsts <- list(T=T)

poisTSInits <- list(lambda0=0,
                     sigma=1,
                     rho=0)

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
set.seed(1)

#call the C++ compiled MCMC model
C_poisTSModelMCMC$run(niter)

#save the MCMC chain (monitored variables) as a matrix
samples <- as.matrix(C_poisTSModelMCMC$mvSamples)
save(samples,file="../Results/MCMC_Chains/mcmc_model_0_whole_Classic.RData")

###
###get log_prob col indeces
print(C_poisTSModelMCMC$calculateWAIC(nburnin=100000))
cols <- grep("log*",colnames(samples))
print(prod(colMeans(exp(samples[100000:200000,cols]))))