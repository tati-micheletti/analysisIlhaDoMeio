# Loading the land species data

loadData <- function(dataLocation){
  DT <- fread(dataLocation)
  DT[, "Date" :=  as.Date(with(DT,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")]
  DT[, "originDate" :=  as.Date(paste(Year,1,1,sep="-"),"%Y-%m-%d")]
  julDate <- rbindlist(lapply(X = 1:NROW(DT), function(rid){
    currRow <- DT[rid,]
    dt <- currRow[["Date"]]
    orig <- currRow[["originDate"]]
    jd <- data.table(JulianDate = julian.Date(x = dt, 
                                              origin = orig))
    return(jd)
  }))
  DT[, "JulianDate" :=  julDate]
  DT[, "timeSinceStartEradication" :=  julian.Date(x = Date, origin = as.Date(min(unique(DT$Date))))]
  return(DT)
}