AICmatches <- function(previousAIC, currentAIC, tolerance){
  # Is current and previous AIC within the tolerance? Return TRUE or FALSE
  minTolerance <- previousAIC - tolerance
  maxTolerance <- previousAIC + tolerance
  AICmatch <- all(currentAIC < maxTolerance, currentAIC > minTolerance)
  return(AICmatch)
}
