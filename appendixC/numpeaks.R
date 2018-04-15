## Function to calculate the number of waves in a prevalence time series from its derivative
numpeaks <- function(dI){
  ## Calculate the number of times the derivative is zero
  diffs <- diff(sign(dI)) # calculate vector of element-wise differences
  zeros <- sum(diffs != 0) # count how many differences are non-zero (each difference that is != 0 is a zero of the derivative)

  ## Calculate the number of waves from the number of zeros
  if (zeros%%2==1){
    # for an odd number of zeros
    numpeaks <- (zeros+1)/2
  } else {
    # for an even number of zeros
    numpeaks <- zeros/2
  }
  
  return(numpeaks)
}
