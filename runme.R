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
Data <- loadData(dataLocation = "data/countData.csv")
weatherData <- loadWeatherData(dataLocation = file.path(wd, "weatherData.csv"))

fullComb <- data.table(expand.grid(Species = c("Trachylepis atlantica",
                                               "Elaenia ridleyana",
                                               "Johngarthia lagostoma",
                                               "Sula dactylatra"),
                                   Island = c("Ilha do Meio", 
                                              "Ilha Rata")))
fullComb[, shortSpIsland := c("mabuiaMeio", 
                              "elaeniaMeio",
                              "crabMeio", 
                              "maskedBoobyMeio",
                              "mabuiaRata", 
                              "elaeniaRata",
                              "crabRata", 
                              "maskedBoobyRata")]
fullComb[, rbstDsgn := c(rep(FALSE, times = 1), 
                         rep(TRUE, times = 1),
                         rep(FALSE, times = 3),
                         rep(TRUE, times = 1),
                         rep(FALSE, times = 2)
)]
# fullComb[, immgrt := c(rep(FALSE, times = 8),
#                        rep(TRUE, times = 1),
#                        rep(FALSE, times = 3),
#                        rep(TRUE, times = 1)
# )]

finalTableName <- file.path("outputs/finalTable.csv")

runParallel <- FALSE


if (!file.exists(finalTableName)){
  if (runParallel){
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
    
  } else {
    t1 <- Sys.time()
    tb <- rbindlist(lapply(
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
  }
  write.csv(tb, finalTableName)
} else {
  tb2 <- fread(finalTableName)
}
 