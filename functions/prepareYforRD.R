prepareYforRD <- function(y, howManySecondary){
  allDates <- unique(colnames(y))
  dt <- unlist(lapply(seq_along(howManySecondary), function(secOcc){
    howManyNow <- howManySecondary[secOcc]
    howManyMax <- max(howManySecondary)
    howManyNA <- howManyMax-howManyNow
    newY <- c(y[colnames(y) == allDates[secOcc]], rep(NA, times = howManyNA))
    return(newY)
  }))
 return(dt)
}
