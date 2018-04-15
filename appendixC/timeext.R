## Set up function to calculate time to extinction
timeext <- function(R0=2,
                    gamma=0.5,
                    rbeta=0.1,
                    betaF=3,
                    gammaF=0.05,
                    N=1e6,
                    rbetaF=0
                    ){
    
# Set initial conditons
I0 <- 1/N
S0 <- 1-I0
S_F0 <- 0

## Calculate beta
beta <- R0*gamma
  
## Solve ODEs numerically once
tmax <- 300 

soln <- bSIRsolve(tmax=tmax,
                  S0=S0, S_F0=S_F0, I0=I0,
                  beta=beta, beta_F=betaF, gamma_F=gammaF,
                  r_beta=rbeta, gamma=gamma, N=N, rbetaF=rbetaF)

## Solve ODEs numerically until the number of infecteds is below one individual

while (soln$I[nrow(soln)]>1/N){
    tmax <- 1.5*tmax
    
    soln <- bSIRsolve(tmax=tmax,
                      S0=S0, S_F0=S_F0, I0=I0,
                      beta=beta, beta_F=betaF, gamma_F=gammaF,
                      r_beta=rbeta, gamma=gamma, N=N, rbetaF=rbetaF)
}

# bSIR data is cut off at the point where soln$I dips below 1/N
# so the maximal output time is the time of extinction
timeext <- max(soln$time)

return(timeext)
}
