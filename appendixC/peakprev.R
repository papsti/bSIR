# Function to calculate the peak prevalence(s)
peakprev <- function(R0=2,
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

## Now we need to find the peaks 
## They will occur where dI/dt=0, but specifically when
## the derivative goes from positive to negative
## i.e. the difference between elements is negative

diffs <- diff(sign(soln$dI)) # calculate vector of element-wise differences
ind <- which(diffs<0) # find the indices for any negative differences
peakdata <- soln[ind,] # extract rows of ODE solution data where a peak occurs

return(peakdata)
}
