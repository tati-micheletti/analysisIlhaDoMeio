getBestModelsTable <- function(fullModTable = NULL){
fileBestModelTable <- "outputs/bestModelTable.csv"
    if (is.null(fullModTable)){
    if (file.exists(fileBestModelTable)){
      bestModTable <- data.table::fread(fileBestModelTable)
    } else {
      fileFullModelTable <- "outputs/fullModelsTable.csv"
      if (file.exists(fileFullModelTable)){
        fullModTable <- data.table::fread(fileFullModelTable)
        bestModTable <- rbindlist(lapply(unique(fullModTable$species), function(sp){
          rbindlist(lapply(unique(fullModTable$island), function(isl){
            fullModTable[species == sp &
                           island == isl &
                           AIC == min(fullModTable[species == sp &
                                                     island == isl, AIC]), ]
          }))
        }))
      } else stop(paste0("Neither bestModelTable.csv or fullModelsTable.csv have",
                         "been created yet. Please create the object 'fullModTable'",
                         "and try again."))
    }
  } else {
    bestModTable <- rbindlist(lapply(unique(fullModTable$species), function(sp){
      rbindlist(lapply(unique(fullModTable$island), function(isl){
        fullModTable[species == sp &
                       island == isl &
                       AIC == min(fullModTable[species == sp &
                                                 island == isl, AIC]), ]
      }))
    }))
  }
if (!file.exists(fileBestModelTable))
  write.csv(bestModTable, file = fileBestModelTable)
return(bestModTable)
}
