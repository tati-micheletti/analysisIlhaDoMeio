# Shortcut wildlife
if(!require("Require")){
  install.packages("Require")
}
library("Require")
Require("googledrive")
Require("reproducible")
Require("data.table")
Require("unmarked")
Require("tictoc")
Require("qs")
# If any packages were installed here, restart your session and re-run the code!
# Here we list all functions and then source all of them:
allFls <- list.files(path = "./functions", full.names = TRUE)
for (fl in allFls){
  message(paste0("Sourcing the function ", basename(fl))) 
  source(file = fl) 
}
Data <- loadESData(dataLocation = "data/IPA.csv")
mabuiaMeio <- modelTerrestrials(DT = Data, island = "Ilha do Meio", species = "Trachylepis atlantica")
mabuiaRata <- modelTerrestrials(DT = Data, island = "Ilha Rata", species = "Trachylepis atlantica")
elaeniaMeio <- modelTerrestrials(DT = Data, island = "Ilha do Meio", species = "Elaenia ridleyana")
elaeniaRata <- modelTerrestrials(DT = Data, island = "Ilha Rata", species = "Elaenia ridleyana")
crabMeio <- modelCrab(DT = Data, island = "Ilha do Meio", species = "Johngarthia lagostoma")
crabRata <- modelCrab(DT = Data, island = "Ilha Rata", species = "Johngarthia lagostoma")
