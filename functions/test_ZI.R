# TEST FOR ZERO INFLATION

test_ZI <- function(dts){
  pois_data <- dts
  lambda_est <- mean(pois_data)
  p0_tilde <- exp(-lambda_est)
  p0_tilde
  n0 <- sum(1*(!(pois_data >0)))
  n <- length(pois_data)
  # number of observtions 'expected' to be zero
  expectedZeros <- n*p0_tilde
  realZeros <- table(pois_data)["0"]
  message(paste0("Real zeros over expected zeros by a factor of ", 
                 realZeros/expectedZeros)) # More realZeros than expected zeros, likely ZI))
  #now lets perform the JVDB score test 
  numerator <- (n0 -n*p0_tilde)^2
  denominator <- n*p0_tilde*(1-p0_tilde) - n*lambda_est*(p0_tilde^2)
  test_stat <- numerator/denominator
  pvalue <- pchisq(test_stat,df=1, ncp=0, lower.tail=FALSE)
  if (pvalue < 0.05){
    message("Data likely zero inflated")
    return(TRUE)
  } else {
    message("Data likely not zero inflated")
    return(FALSE)
  }
}
