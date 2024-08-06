joinCols <- function(y, DT, numPrimary){
  interval <- NCOL(y)/numPrimary # How many secondary periods we have
  firstColsToJoin <- seq(from = 1, to = NCOL(y), by = interval)
  DTnew <- data.table(DT)
  newDT <- lapply(firstColsToJoin, function(i){
    colsToJoin <- colnames(DT)[i:(i+1)]
    DT2 <- DTnew[, .SD, .SDcols = colsToJoin]
    DT2[, colnames(DT2)[1] := rowMeans(.SD), .SDcols = names(DT2)]
    DT2[,names(DT2)[NCOL(DT2)] := NULL]
    return(DT2)
  })
  newDT <- setDT(unlist(newDT, recursive = FALSE),check.names = TRUE)[]
  return(newDT)
}

# Alternative trial for robust design informing how many secondary periods
# However, these are not always the same amount. So that's a problem unless I 
# fill it up with NA's.
# joinCols <- function(DT, numPrimary){
#   allDates <- colnames(DT)
#   newDT <- lapply(allDates, function(i){
#     DT2 <- data.table(DT[, i])
#     names(DT2) <- i
#     DT2[, colnames(DT2)[1] := rowMeans(.SD), .SDcols = names(DT2)]
#     return(DT2)
#   })
#   newDT <- setDT(unlist(newDT, recursive = FALSE),check.names = TRUE)[]
#   return(newDT)
# }