# TODO
# [X] Check all models for sense. Double check the estimates produced by the last 2 K models. 
# Gamma doesn't change for any of the laste 2 K's, even though  AIC is not stable (<0.8 difference). Will go with this! 
# [X] Re-run needed models: Elaenia Meio and Rata
# [X] Re-ran Rata without robust design. Growth of 1.05 makes sense... Re-running currently Meio with Robust design.
# Meio with Robust design still shows negative negative correlation with TSE... sigh. Not sure what to do here...
# Rata without robust design results saved .
# [X] Without robust design elaenia not better: try making robust design again. -- not good...
  # NOTE -- Without Robust design the gamma for Elaenia on Meio is wacked with Immigration. 
# [X] We will try a new variable, using amount of rain as opposed to Julian Date. It might be better for all models. 
# [X] MAKE SURE TO ADD TO THE NEW TABLE'S NAMES THE PRESENCE OF IMMIGRATION AND RD SOMEHOW?
# [ ] To be comparable, the models need to be the same relative to Robust Design and Immigration!
# [ ] Redo the final tables!! (Should work just re-running the code as all the models have been built already and saved...)

library("qs")
library("unmarked")
# Elaenia
emMod <- qread("~/projects/analysisIlhaDoMeio/outputs/elaenia_Meio_bestModel.qs")
emP <- predict(obj = emMod[[length(emMod)]], type = "gamma")
erMod <- qread("~/projects/analysisIlhaDoMeio/outputs/elaenia_Rata_bestModel.qs")
erP <- predict(obj = erMod[[length(erMod)]], type = "gamma")


# Crab
cmMod <- qread("~/projects/analysisIlhaDoMeio/outputs/crab_Meio_bestModel.qs")
cmP <- predict(obj = cmMod[[length(cmMod)]], type = "gamma")
crMod <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/crab_Rata_K1650.qs")
crP <- predict(obj = crMod, type = "gamma")
crMod2 <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/crab_Rata_K1450.qs")
crP2 <- predict(obj = crMod, type = "gamma")
crPred <- data.table::data.table(K1650 = crP[,"Predicted"],
                                 K1450 = crP2[,"Predicted"])

# mabuia
mmMod <- qread("~/projects/analysisIlhaDoMeio/outputs/mabuia_Meio_bestModel.qs")
mmP <- predict(obj = mmMod[[length(mmMod)]], type = "gamma")
mrMod <- qread("~/projects/analysisIlhaDoMeio/outputs/mabuia_Rata_bestModel.qs")
mrP <- predict(obj = mrMod[[length(mrMod)]], type = "gamma")

# booby
bmMod <- qread("~/projects/analysisIlhaDoMeio/outputs/maskedBooby_Meio_bestModel.qs")
bmP <- predict(obj = bmMod[[length(bmMod)]], type = "gamma")
brMod <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/maskedBooby_Rata_K1852.qs")
brP <- predict(obj = brMod, type = "gamma")
brMod2 <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/maskedBooby_Rata_K1652.qs")
brP2 <- predict(obj = brMod, type = "gamma")
brPred <- data.table::data.table(K1852 = brP[,"Predicted"],
                                 K1652 = brP2[,"Predicted"])

#1250, 1050
library("qs")
library("unmarked")
emMod <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/elaenia_Meio_K1250.qs")
emP <- predict(obj = emMod, type = "gamma")
emMod2 <- qread("~/projects/analysisIlhaDoMeio/outputs/individualModels/elaenia_Meio_K1050.qs")
emP2 <- predict(obj = emMod2, type = "gamma")
emPred <- data.table::data.table(K1250 = emP[,"Predicted"],
                                 K1050 = emP2[,"Predicted"])


emBMod <- qread("~/projects/analysisIlhaDoMeio/outputs/elaenia_Meio_bestModel.qs")
emBPred <- predict(obj = emBMod$K_250, type = "gamma")

erBMod <- qread("~/projects/analysisIlhaDoMeio/outputs/elaenia_Rata_bestModel.qs")
erBPred <- predict(obj = erMod$K_1850, type = "gamma")

eMiPred <- predict(obj = eMi, type = "gamma")
