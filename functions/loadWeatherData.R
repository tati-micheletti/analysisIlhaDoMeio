loadWeatherData <- function(dataLocation){
  DT <- fread(dataLocation)
  names(DT)[names(DT) == "datetime"] <- "Date"
  setDT(DT)[, c("Year", "Month", "Day") := tstrsplit(Date, "-")]
  DT[, "Date" :=  as.Date(Date,"%Y-%m-%d")]
  return(DT)
  
}
