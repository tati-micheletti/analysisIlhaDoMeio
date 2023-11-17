joinCols <- function(y, DT, numPrimary){
  interval <- NCOL(y)/numPrimary # How many secondary periods we have
  firstColsToJoin <- seq(from = 1, to = NCOL(y), by = interval)
  newDT <- lapply(firstColsToJoin, function(i){
    colsToJoin <- colnames(DT)[i:(i+1)]
    DT2 <- data.table(DT[, colsToJoin])
    DT2[, colnames(DT2)[1] := rowMeans(.SD), .SDcols = names(DT2)]
    DT2[,names(DT2)[NCOL(DT2)] := NULL]
    return(DT2)
  })
  newDT <- setDT(unlist(newDT, recursive = FALSE),check.names = TRUE)[]
  return(newDT)
}
