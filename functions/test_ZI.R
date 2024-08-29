# TEST FOR ZERO INFLATION

test_ZI <- function(dts, percTolerance = 0.20){
  pois_data <- dts
  lambda_est <- mean(pois_data)
  p0_tilde <- exp(-lambda_est)
  n0 <- sum(1*(!(pois_data >0)))
  lng <- length(pois_data)
  # number of observtions 'expected' to be zero
  expectedZeros <- lng*p0_tilde
  realZeros <- table(pois_data)["0"]
  message(paste0("Real zeros over expected zeros by a factor of ", 
                 round(realZeros/expectedZeros, digits = 2))) # More realZeros than expected zeros, likely ZI))
  #now lets perform the JVDB score test 
  numerator <- (n0 -lng*p0_tilde)^2
  denominator <- lng*p0_tilde*(1-p0_tilde) - lng*lambda_est*(p0_tilde^2)
  test_stat <- numerator/denominator
  pvalue <- pchisq(test_stat,df=1, ncp=0, lower.tail=FALSE)
  perc <- realZeros/lng
  message(paste0("Percentage of zeros in dataset  is ", 
                 round(perc, digits = 2), " tolerance is ", 
                 round(percTolerance, digits = 2))) # High percentage of zeros (rule of thumb 20-50%), likely ZI))
  if (all(perc > percTolerance, pvalue < 0.05)){
    message("Data likely zero inflated")
    return(TRUE)
  } else {
    message("Data likely not zero inflated")
    return(FALSE)
  }
}
