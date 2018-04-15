## Set up parms.df (parameter dataframe)
R0.vec <- c(rep(2,7)) # basic reproduction number
beta_F.vec <- c(0.25, 0, 1, 2.5, 5, 3, 3) # transmission of fear
gamma_F.vec <- c(rep(0.5,5), 0.1, 0.05) # recovery rate from fear
r_beta.vec <- c(rep(0.5,5), 0.1, 0.1) # constant of proportionality 
# between transmission rate for susceptibles and fearful susceptibles
gamma.vec <- c(rep(0.1,6), 0.5) # recovery rate from the disease
N.vec <- c(rep(1e6,7)) # population size
# transmission rate of the disease from the S compartment
# and the S_F compartment

## Define transmission term
beta.vec <- R0.vec*gamma.vec

## Set initial conditions (state variables are in numbers of individuals)
I0.vec <- rep(1, length(R0.vec))
S0.vec <- N.vec-I0.vec
S_F0.vec <- rep(0, length(R0.vec)) # no one is initially fearful

## Set max integration time
tmax <- 300

## Set up parameter array
parm.df <- data.frame(R0=R0.vec,
                      beta=beta.vec,
                      beta_F=beta_F.vec,
                      gamma_F=gamma_F.vec,
                      r_beta=r_beta.vec,
                      gamma=gamma.vec,
                      N=N.vec,
                      S0=S0.vec,
                      S_F0=S_F0.vec,
                      I0=I0.vec)
for(i in 1:7){
  parmsrow <- i
  
  ## Set parameters
  R0 <- parm.df$R0[parmsrow]
  beta <- parm.df$beta[parmsrow]
  beta_F <- parm.df$beta_F[parmsrow]
  gamma_F <- parm.df$gamma_F[parmsrow]
  r_beta <- parm.df$r_beta[parmsrow]
  gamma <- parm.df$gamma[parmsrow]
  N <- parm.df$N[parmsrow]
  S0 <- parm.df$S0[parmsrow]
  S_F0 <- parm.df$S_F0[parmsrow]
  I0 <- parm.df$I0[parmsrow]
  
  ## Generate data
  bSIRdata <- bSIRsolve(tmax=tmax,
                        S0=S0, S_F0=S_F0, I0=I0,
                        beta=beta, beta_F=beta_F, gamma_F=gamma_F,
                        r_beta=r_beta, gamma=gamma)
  
  ## Relabel the data with a fitsir-friendly format
  bSIRdata <- data.frame(tvec=bSIRdata$time, count=bSIRdata$I)
  
  ## Store parameter set used for these data
  parms <- parm.df[parmsrow,]
  
  ## Save objects to .Rdata file
  save(list=c("bSIRdata", "parms"), 
       file=paste0("data/bSIR-parmsrow", parmsrow, ".Rdata"))
}