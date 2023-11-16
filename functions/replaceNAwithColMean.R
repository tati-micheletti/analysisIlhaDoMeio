replaceNAwithColMean  <- function(mx){
  dt <- as.data.table(mx)
  nm <- names(dt)[colSums(is.na(dt)) != 0]
  for(clnm in nm){
    set(dt, i = which(is.na(dt[[clnm]])), j = clnm, 
        value = na.omit(unique(dt[, get(clnm)], na.rm = TRUE)))
  }
  return(as.matrix(dt))
}
