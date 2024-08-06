loadWeatherData <- function(dataLocation){
  DT <- fread(dataLocation)
  names(DT)[names(DT) == "datetime"] <- "Date"
  setDT(DT)[, c("Month", "Day", "Year") := tstrsplit(Date, "/")]
  DT[, "Date" :=  as.Date(Date,"%m/%d/%Y")]
  return(DT)
}
