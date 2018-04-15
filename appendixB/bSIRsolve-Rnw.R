## Load package to numerically solve ODEs
library("deSolve")

## Define a function to numerically solve the bSIR equations
## with a C gradient for speed
bSIRsolve <- function(tmax, # max integration time, in days
                      S0, S_F0, I0, # initial conditions in proportions of the population
                      beta, beta_F, gamma_F, r_beta, gamma, N, # parameters
                      rbetaF=0 # rbetaF value for prev- and fear-based spread (if 0, then there is just prev-based fear)
                      ){
  
  ## Set integration parameters
  stepsize <- 1/24
  
  # Solve the odes numerically using the previously defined R gradient
  Z0 <- I0 # set initial cumulative case count equal to the initial number of infecteds

    ## Load DLL
    dyn.load("~/Dropbox/bSIR/scripts/bSIR_PF.so") 
    ## Solve ODE
    soln <- ode(y=c(S=S0,S_F=S_F0,I=I0,Z=Z0),
                times=seq(0, tmax, by=stepsize),
                func="derivs",
                parms=c(beta=beta, beta_F=beta_F, gamma_F=gamma_F, r_beta=r_beta, gamma=gamma, rbetaF=rbetaF),
                dllname = "bSIR_PF",
                initfunc = "initmod",
                nout = 2, outnames = c("R", "dI"))
    soln <- as.data.frame(soln)
  
  ## If prevalence falls below one individual
  ## stop returning data
  ind <- which(soln$I<1/N)
  if (length(ind)>0){
  soln <- soln[1:min(ind),]
  }

  return(soln)
}
