# Check individual models
library("data.table")

# allSp <- c("Trachylepis atlantica", 
#            "Elaenia ridleyana", 
#            "Johngarthia lagostoma", 
#            "Sula dactylatra")
# allIsl <- c("Ilha Rata", "Ilha do Meio")
# 
# mM <- tb2[island == "Ilha do Meio" &
#             species ==  "Trachylepis atlantica",]
# 
# sp <- c("crab", "mabuia", "maskedBooby", "elaenia")
# isl <- c("Meio", "Rata")
# tp <- c("individualModels_", 
#           "individualModels_Immig", 
#           "individualModels_Immig_RbstDsgn", 
#           "individualModels_RbstDsgn")

best_mod <- function(species, island, modType){
outDir <- file.path(getwd(), "outputs")
  allMds <- list.files(path = outDir, 
            pattern = "bestModel", 
            recursive = TRUE)
mod <- qs::qread(file.path(outDir, modType, paste0(species, "_", island,"_bestModel.qs")))
}

DT <-  tb2 <- fread(file.path("outputs/finalTable.csv"))

# Mabuia
mM <- best_mod(species = "mabuia", island = "Meio", modType = "individualModels_")
pmM <- unmarked::predict(mM[[length(mM)]], type = "gamma")
mR <- best_mod(species = "mabuia", island = "Rata", modType = "individualModels_")
pmR <- unmarked::predict(mR[[length(mR)]], type = "gamma")

# Crab
cM <- best_mod(species = "crab", island = "Meio", modType = "individualModels_")
pcM <- unmarked::predict(cM[[length(cM)]], type = "gamma")
cR <- best_mod(species = "crab", island = "Rata", modType = "individualModels_")
pcR <- unmarked::predict(cR[[length(cR)]], type = "gamma")

# Masked Booby 
mbM <- best_mod(species = "maskedBooby", island = "Meio", modType = "individualModels_")
pmbM <- unmarked::predict(mbM[[length(mbM)]], type = "gamma")
DT[island == "Ilha do Meio" & species == "Sula dactylatra" & 
     AIC == min(DT[island == "Ilha do Meio" & species == "Sula dactylatra", AIC]), ]

mbR <- best_mod(species = "maskedBooby", island = "Rata", modType = "individualModels_")
pmbR <-unmarked::predict(mbR[[length(mbR)]], type = "gamma")
DT[island == "Ilha Rata" & species == "Sula dactylatra" & 
     AIC == min(DT[island == "Ilha Rata" & species == "Sula dactylatra", AIC]), ]

# Masked Booby Immigration
mbMI <- best_mod(species = "maskedBooby", island = "Meio", modType = "individualModels_Immig")
pmbMI <- unmarked::predict(mbMI[[length(mbMI)]], type = "gamma")
DT[island == "Ilha do Meio" & species == "Sula dactylatra" & 
     AIC == min(DT[island == "Ilha do Meio" & species == "Sula dactylatra", AIC]), ]

mbRI <- best_mod(species = "maskedBooby", island = "Rata", modType = "individualModels_Immig")
pmbRI <- unmarked::predict(mbRI[[length(mbRI)]], type = "gamma")
DT[island == "Ilha Rata" & species == "Sula dactylatra" & 
     AIC == min(DT[island == "Ilha Rata" & species == "Sula dactylatra", AIC]), ]

# Elaenia
eM <- best_mod(species = "elaenia", island = "Meio", modType = "individualModels_RbstDsgn")
peM <- unmarked::predict(eM[[length(eM)]], type = "gamma")
DT[island == "Ilha do Meio" & species == "Elaenia ridleyana" & 
     AIC == min(DT[island == "Ilha do Meio" & species == "Elaenia ridleyana", AIC]), ]

eR <- best_mod(species = "elaenia", island = "Rata", modType = "individualModels_RbstDsgn")
peR <- unmarked::predict(eR[[length(eR)]], type = "gamma")
DT[island == "Ilha Rata" & species == "Elaenia ridleyana" & 
     AIC == min(DT[island == "Ilha Rata" & species == "Elaenia ridleyana", AIC]), ]

# Elaenia Immigration
eMI <- best_mod(species = "elaenia", island = "Meio", modType = "individualModels_Immig_RbstDsgn")
peMI <- unmarked::predict(eM[[length(eM)]], type = "gamma")
DT[island == "Ilha do Meio" & species == "Elaenia ridleyana" & 
     AIC == min(DT[island == "Ilha do Meio" & species == "Elaenia ridleyana", AIC]), ]

eRI <- best_mod(species = "elaenia", island = "Rata", modType = "individualModels_Immig_RbstDsgn")
peRI <- unmarked::predict(eR[[length(eR)]], type = "gamma")
DT[island == "Ilha Rata" & species == "Elaenia ridleyana" & 
     AIC == min(DT[island == "Ilha Rata" & species == "Elaenia ridleyana", AIC]), ]

### PLOTS

makePlot <- function(predMeio, predRata, Species){
  eDt <- data.table(x = 1:length(predMeio$Predicted), 
                    y = predMeio$Predicted-predRata$Predicted)
  plot(main = Species, x = 1:length(predMeio$Predicted), y = predMeio$Predicted-predRata$Predicted, 
       xlab = "Time", ylab = "Growth (lambda) Difference")
  eFit <- lm(y ~ x, data = eDt)
  abline(eFit)
  abline(h = 0, col = "red")
}

makePlot(predMeio = peM, predRata = peR, Species = "Elaenia") # Elaenia: 12 occasions
makePlot(predMeio = pmbMI[1:3,], predRata = pmbRI[1:3,], Species = "Masked Booby") # Booby: 3 occasions
makePlot(predMeio = pmM[1:4,], predRata = pmR[1:4,], Species = "Mabuia") # Mabuia: 4 occasions
makePlot(predMeio = pcM[1:10,], predRata = pcR[1:10,], Species = "Crab") # Mabuia: 4 occasions

