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
wd <- file.path("./data")
# If any packages were installed here, restart your session and re-run the code!
# Here we list all functions and then source all of them:
allFls <- list.files(path = "./functions", full.names = TRUE)
for (fl in allFls){
  message(paste0("Sourcing the function ", basename(fl))) 
  source(file = fl) 
}
# Data <- loadData(dataLocation = "data/dataCOUNTS.csv")
# mabuiaMeio <- modelWildlife(DT = Data, island = "Ilha do Meio", species = "Trachylepis atlantica")
# mabuiaRata <- modelWildlife(DT = Data, island = "Ilha Rata", species = "Trachylepis atlantica")
# elaeniaMeio <- modelWildlife(DT = Data, island = "Ilha do Meio", species = "Elaenia ridleyana")
# elaeniaRata <- modelWildlife(DT = Data, island = "Ilha Rata", species = "Elaenia ridleyana")
changeInPopulation <- calculatePopulationChange(species = c(
  "Trachylepis atlantica",
  "Elaenia ridleyana"))

blackRatData <- invisible(getRatDensity())
plotGrowth <- plotPopulationGrowth(wildlifeDataset = changeInPopulation, 
                                   ratsDataset = blackRatData)
# crabMeio <- modelWildlife(DT = Data, island = "Ilha do Meio", species = "Johngarthia lagostoma")
# crabRata <- modelWildlife(DT = Data, island = "Ilha Rata", species = "Johngarthia lagostoma")
# boobyMeio <- modelWildlife(DT = Data, island = "Ilha do Meio", species = "Sula dactylatra")
# boobyRata <- modelWildlife(DT = Data, island = "Ilha Rata", species = "Sula dactylatra")

