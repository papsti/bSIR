fitsirplot <- function(paramset=2,
                       bSIRlwd=12,
                       plotcex=0.75,
                       margins=c(2.5,2.5,0.1,0.1)+0.1
){
  
  ## Load pre-generated bSIR data
  load(paste0("data/bSIR-parmsrow", paramset, ".Rdata"))
  
  ## Source files
  library("fitsir")
  library("bbmle")
  
  ## Rename parameters for cleaner code
  R0 <- parms$R0
  beta <- parms$beta
  beta_F <- parms$beta_F
  gamma_F <- parms$gamma_F
  r_beta <- parms$r_beta
  gamma <- parms$gamma
  N <- parms$N
  S0 <- parms$S0
  S_F0 <- parms$S_F0
  I0 <- parms$I0
  
  ## Fit calculation
  m1 <- fitsir(data=bSIRdata, start=startfun(log.beta=log(beta),
                                             log.gamma=log(gamma),
                                             log.N=log(N),
                                             logit.i=qlogis(I0/N)))
  
  ## Simulate SIR with fitted parameters
  SIRfitdata <- with(bSIRdata, SIR.detsim(tvec,trans.pars(coef(m1))))
  
  ## Set plotting parameters
  bSIRcol <- "#FF8000"
  SIRcol <- "blue"
  SIRlwd <- bSIRlwd/3
  
  ## Set up plotting region
  #par(mar=c(5,5,1,8.25)+0.1)
  par(mar=margins)
  
  ## Set mult factor to expand plot limits
  limmult <- 1.05
  
  ## Set mult factor to position legend in margin
  xlegmult <- 1.05
  ylegmult <- 1.13
  
  ## Plot bSIR "data"
  plot((count/N)~tvec, data=bSIRdata,
       type="l", xaxs="i", 
       xlab="", ylab="",
       col=bSIRcol, lwd=bSIRlwd,
       ylim=c(0, max(count/N)),
       cex.axis=plotcex, cex.lab=plotcex, las=1)
  
  ## Add fitted SIR solution
  lines(bSIRdata$tvec, SIRfitdata/N, col=SIRcol, lwd=SIRlwd)
  
  ## Log-Likelihood
  loglikelihood <- sprintf("$%g$", as.numeric(logLik(m1)))
  legend(0.6*max(bSIRdata$tvec), 1.05*max(bSIRdata$count/N),
         #xlegmult*max(bSIRdata$tvec), -0.03*max(bSIRdata$count),
         legend=loglikelihood,
         bty="n", title="\\textbf{Log-Likelihood}", xpd=TRUE,
         cex=plotcex, title.adj=0)
  
  ## Parameter set label
  parsetlab <- paste0("Parameter Set ", paramset-1)
  title(main=paste0("$\\textbf{", parsetlab, "}$"), line=1)
}