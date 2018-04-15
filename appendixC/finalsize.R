## Set up function to calculate final size
finalsize <- function(R0=2,
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
  
## Solve odes numerically initially
tmax <- 300
soln <- bSIRsolve(tmax=tmax,
                  S0=S0, S_F0=S_F0, I0=I0,
                  beta=beta, beta_F=betaF, gamma_F=gammaF,
                  r_beta=rbeta, gamma=gamma, N=N, rbetaF=rbetaF)

## Continue integrating until the epidemic burns out
while (soln$I[nrow(soln)]>1/N){
    tmax <- 1.5*tmax
    
    soln <- bSIRsolve(tmax=tmax,
                  S0=S0, S_F0=S_F0, I0=I0,
                  beta=beta, beta_F=betaF, gamma_F=gammaF,
                  r_beta=rbeta, gamma=gamma, N=N, rbetaF=rbetaF)
    }

## bSIRsolve output is cut off when I<1/N, so the last (maximal) value of Z is the final size
finalsize <- max(soln$Z)

return(finalsize)
}
