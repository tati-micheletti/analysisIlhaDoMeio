# Shortcut wildlife
if(!require("Require")){
  install.packages("Require")
}
library("Require")
Require("googledrive")
Require("reproducible")
Require("data.table")
Require("unmarked")
Require("terra")
Require("secr")
Require("tictoc")
Require("qs")
Require("ggplot2")
Require("future")
Require("future.apply")

wd <- file.path("./data")
# If any packages were installed here, restart your session and re-run the code!
# Here we list all functions and then source all of them:
allFls <- list.files(path = "./functions", full.names = TRUE)
for (fl in allFls){
  message(paste0("Sourcing the function ", basename(fl))) 
  source(file = fl) 
}
Data <- loadData(dataLocation = "data/dataCOUNTS.csv")
weatherData <- loadWeatherData(dataLocation = file.path(wd, "weatherNoronhaMar17Apr18.csv"))

finalTableName <- file.path("data/finalTable.csv")
if (!file.exists(finalTableName)){
  fullComb <- data.table(expand.grid(Species = c("Trachylepis atlantica",
                                      "Elaenia ridleyana",
                                      "Elaenia ridleyana",
                                      "Elaenia ridleyana",
                                      "Elaenia ridleyana",
                                      "Johngarthia lagostoma",
                                      "Sula dactylatra",
                                      "Sula dactylatra"),
                          Island = c("Ilha do Meio", 
                                     "Ilha Rata")))
  fullComb[, shortSpIsland := c("mabuiaMeio", 
                                "elaeniaMeio",
                                "elaeniaMeio",
                                "elaeniaMeio",
                                "elaeniaMeio",
                                "crabMeio", 
                                "maskedBoobyMeio",
                                "maskedBoobyMeio",
                                "mabuiaRata", 
                                "elaeniaRata",
                                "elaeniaRata",
                                "elaeniaRata",
                                "elaeniaRata",
                                "crabRata", 
                                "maskedBoobyRata",
                                "maskedBoobyRata")]
  fullComb[, rbstDsgn := c(rep(FALSE, times = 3), 
                           rep(TRUE, times = 2),
                           rep(FALSE, times = 4),
                           rep(TRUE, times = 2),
                           rep(FALSE, times = 5)
                           )]
  fullComb[, immgrt := c(rep(FALSE, times = 2),
                           rep(TRUE, times = 2),
                           rep(FALSE, times = 3),
                         TRUE,
                         FALSE,TRUE,FALSE,TRUE, # Rata
                           rep(FALSE, times = 3),
                           rep(TRUE, times = 1)
  )]
  
  # TRY TO RUN IN PARALLEL! Check how we used future in birds stuff
  nCoresNeeded <- NROW(fullComb)
  nCoresAvail <- min(parallel::detectCores() - 3, getOption("NCONNECTIONS", 120L))
  nCores2Use <- min(nCoresNeeded, nCoresAvail)
  plan("multicore", workers = nCores2Use)
  t1 <- Sys.time()
  tb <- rbindlist(future_lapply(
    1:NROW(fullComb),
    function(index) {
      modelAndSummary(
        ROW = index,
        fullComb = fullComb,
        Data = Data,
        weatherData = weatherData
      )
    }
  ))
  plan("sequential")
  write.csv(tb, finalTableName)
} else {
  tb2 <- fread(finalTableName)
}
  
# blackRatData <- invisible(getRatDensity())
# plotGrowth <- plotPopulationGrowth(wildlifeDataset = changeInPopulation, 
#                                    ratsDataset = blackRatData)

