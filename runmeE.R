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
Data <- loadData(dataLocation = "data/finalDataCOUNTS.csv")
weatherData <- loadWeatherData(dataLocation = file.path(wd, "weatherNoronhaMar17Nov23.csv"))

fullComb <- data.table(expand.grid(Species = c("Trachylepis atlantica",
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
                              "crabMeio", 
                              "maskedBoobyMeio",
                              "maskedBoobyMeio",
                              "mabuiaRata", 
                              "elaeniaRata",
                              "elaeniaRata",
                              "crabRata", 
                              "maskedBoobyRata",
                              "maskedBoobyRata")]
fullComb[, rbstDsgn := c(rep(FALSE, times = 1), 
                         rep(TRUE, times = 2),
                         rep(FALSE, times = 4),
                         rep(TRUE, times = 2),
                         rep(FALSE, times = 3)
)]
fullComb[, immgrt := c(rep(FALSE, times = 2),# Meio
                       rep(TRUE, times = 1),
                       rep(FALSE, times = 2),
                       rep(TRUE, times = 1),
                       rep(FALSE, times = 2),# Rata
                       rep(TRUE, times = 1),
                       rep(FALSE, times = 2),
                       rep(TRUE, times = 1)
)]

fullComb2 <- fullComb[Species == "Elaenia ridleyana",]

nCoresNeeded <- NROW(fullComb2)
nCoresAvail <- min(parallel::detectCores() - 3, getOption("NCONNECTIONS", 120L))
nCores2Use <- min(nCoresNeeded, nCoresAvail)
plan("multicore", workers = nCores2Use)
t1 <- Sys.time()
tb <- rbindlist(future_lapply(
  1:NROW(fullComb2),
  function(index) {
    modelAndSummary(
      ROW = index,
      fullComb = fullComb2,
      Data = Data,
      weatherData = weatherData,
      rerunModels = TRUE,
      reRunK = TRUE,
      rewriteTable = TRUE
    )
  }
))
plan("sequential")
