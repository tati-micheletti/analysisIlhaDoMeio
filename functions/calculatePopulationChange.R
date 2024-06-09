calculatePopulationChange <- function(species, island,
                                      envir = environment(),
                                      rewriteTable = FALSE,
                                      RobustDesign,
                                      shortName){
  if (any(length(species)>1,
          length(island)>1,
          length(shortName)>1)) stop("Please provide only one species and one island")
  outputsTable <- file.path("outputs", paste0("outputsTable_",shortName,".csv"))
  if (any(!file.exists(outputsTable),
          rewriteTable)) {  
    listSpecies <- rbindlist(lapply(species, function(sp){
      listIsland <- rbindlist(lapply(island, function(isla){
        spShort <- unlist(sapply(sp, switch,
                                 "Elaenia ridleyana" = "elaenia",
                                 "Trachylepis atlantica" = "mabuia",
                                 "Johngarthia lagostoma" = "crab",
                                 "Sula dactylatra" = "maskedBooby", 
                                 USE.NAMES = FALSE))
        islandShort <- unlist(sapply(isla, switch,
                                           "Ilha Rata" = "Rata",
                                           "Ilha do Meio" = "Meio", 
                                     USE.NAMES = FALSE))
        if (islandShort == "Meio"){
          cyclesEach <- as.numeric(sapply(spShort, switch,
                                          "elaenia" = if (RobustDesign) c(1,2,4,5) else c(1,1.1,2,2.1,4,4.1,5,5.1), 
                                          "mabuia" = c(1,2,2.1,4,5), #c(1,2,4,4.1,5),
                                          "crab" = c(1,1.1,2,2.1,4,4.1,5),#c(1,2,2.1,4,4.1,5,5.1,6),
                                          "maskedBooby" = c(1,2,4,5)))
        } else {
          cyclesEach <- as.numeric(sapply(spShort, switch,
                                          "elaenia" = if (RobustDesign) c(1,2,4,5) else c(1,1.1,2,2.1,4,4.1,5,5.1),
                                          "mabuia" = c(1,2,2.1,4,5,5.1), # c(1,2,4,4.1,5,5.1),
                                          "crab" = c(1,2,2.1,2.2,4,5),# c(1,2,4,4.1,4.2,5,6),
                                          "maskedBooby" = c(1,2,4,5)))
        }
        resIsla <- get(paste0(spShort, islandShort), envir = envir)
        listParam <- rbindlist(lapply(c("populationGrowthRate", 
                                        "detectionProbability",
                                        "immigration",
                                        "initialPopulationDensity"), function(param){
                                          message(paste0("Creating data for ", param, " for ", spShort, " for ", isla))
                                          tr <- resIsla$samplingOccasions*resIsla$numberOfSites
                                          # POPULATION GROWTH
                                          if (param %in% c("populationGrowthRate", "initialPopulationDensity")){
                                            lng <- length(unique(resIsla$populationGrowthRate[,"Predicted"]))
                                            if (!any(lng == (resIsla$samplingOccasions-1), lng == 1)){
                                              message("Sites differ for population growth rate. Check.")
                                              browser()
                                            }
                                            popGrRate <- transpose(resIsla$populationGrowthRate[1:as.numeric(resIsla$samplingOccasions-1),])
                                            popGrRate <- cbind(popGrRate, rep(NA, NROW(popGrRate)))
                                            # colnames(popGrRate) <- c(paste0("Cycle_", 1:resIsla$samplingOccasions))
                                            colnames(popGrRate) <- c(paste0("Cycle_", round(cyclesEach, 1)))
                                            rownames(popGrRate) <- c("Predicted", "SE", "lower", "upper")
                                            estimate <- rep(as.numeric(popGrRate["Predicted",]), each = resIsla$numberOfSites)
                                            SE <- rep(as.numeric(popGrRate["SE",]), each = resIsla$numberOfSites)
                                            LCI <- rep(as.numeric(popGrRate["lower",]), each = resIsla$numberOfSites)
                                            UCI <- rep(as.numeric(popGrRate["upper",]), each = resIsla$numberOfSites)
                                          }
                                          # DETECTION PROBABILITY
                                          if (param == "detectionProbability"){
                                            # Need to reorganize the data because it comes as row lead (i.e., cycle 1 to 5 for site 1, then cycle 1 to 5 for site 2 and so on)
                                            estimate <- as.numeric(matrix(resIsla$detectionProbability[,"Predicted"], nrow = resIsla$numberOfSites, byrow = TRUE))
                                            SE <-  as.numeric(matrix(resIsla$detectionProbability[,"SE"], nrow = resIsla$numberOfSites, byrow = TRUE))
                                            LCI <- as.numeric(matrix(resIsla$detectionProbability[,"lower"], nrow = resIsla$numberOfSites, byrow = TRUE))
                                            UCI <- as.numeric(matrix(resIsla$detectionProbability[,"upper"], nrow = resIsla$numberOfSites, byrow = TRUE))
                                          }
                                          # IMMIGRATION
                                          if (param == "immigration"){
                                            if (length(resIsla$immigration) == 1){
                                              if (!is.na(resIsla$immigration)) stop("Length immigration is 1 but not NA. Debug.")
                                              estimate <- SE <- UCI <- LCI <- rep(NA, times = tr)
                                            } else {
                                              imm <- transpose(resIsla$immigration[1:as.numeric(resIsla$samplingOccasions-1),])
                                              imm <- cbind(rep(NA, NROW(imm)), imm)
                                              colnames(imm) <- c(paste0("Cycle_", 1:resIsla$samplingOccasions))
                                              rownames(imm) <- c("Predicted", "SE", "lower", "upper")
                                              estimate <- rep(as.numeric(imm["Predicted",]), each = resIsla$numberOfSites)
                                              SE <- rep(as.numeric(imm["SE",]), each = resIsla$numberOfSites)
                                              LCI <- rep(as.numeric(imm["lower",]), each = resIsla$numberOfSites)
                                              UCI <- rep(as.numeric(imm["upper",]), each = resIsla$numberOfSites)
                                            } 
                                          }
                                          # AVERAGE DENSITY
                                          Cover <- as.character(resIsla$cover)
                                          Order <- 1:resIsla$numberOfSites
                                          dtCO <- data.table(cover = Cover,
                                                             order = Order)
                                          pop <- merge(resIsla$initialPopulationDensity[, c("cover", "densityM")],
                                                       dtCO, by = "cover")
                                          setkey(pop, "order")
                                          if (param == "initialPopulationDensity"){
                                            newD <- estimate <- pop$densityM
                                            newSE <- SE <- pop$densityM # Even though we have SE, lower and upper in  
                                            newLCI <- LCI <- pop$densityM # initial pop density, we can't use it for calculating density.
                                            newUCI <- UCI <- pop$densityM
                                            for (occ in 1:(length(colnames(popGrRate))-1)){
                                            # for (occ in 1:(resIsla$samplingOccasions-1)) {
                                            # for (occ in 1:length(colnames(popGrRate))) {
                                              # Density
                                              # newD <- newD*popGrRate["Predicted",paste0("Cycle_", occ)]
                                              newD <- newD*popGrRate["Predicted",colnames(popGrRate)[occ]] # Make sure to check this is not constant through time...
                                              estimate <- c(estimate, newD)
                                              # SE
                                              newSE <- newSE*popGrRate["SE",colnames(popGrRate)[occ]] # Make sure to check this is not constant through time...
                                              SE <- c(SE, newSE)
                                              # LCI
                                              newLCI <- newLCI*popGrRate["lower",colnames(popGrRate)[occ]] # Make sure to check this is not constant through time...
                                              LCI <- c(LCI, newLCI)
                                              # UCI
                                              newUCI <- newUCI*popGrRate["upper",colnames(popGrRate)[occ]] # Make sure to check this is not constant through time...
                                              UCI <- c(UCI, newUCI)
                                            }
                                            # Check that the length of newD matches the expected length
                                            if (length(estimate) != tr){
                                              message(paste0("Length of PopulationDensity doesn't match expected length: ",
                                                          tr,". Debug."))
                                              browser()
                                            }
                                          }
                                          # WE HAVE THE 40 SITES then the 5 CYCLES:
                                          # cycle 1: sites 1 to 40,
                                          # cycle 2: sites 1 to 40...
                                          res <- data.table(species = rep(resIsla$species, times = tr),
                                                            island = rep(isla, times = tr),
                                                            cycle = rep(cyclesEach, each = resIsla$numberOfSites),
                                                            site = rep(1:resIsla$numberOfSites, times = resIsla$samplingOccasions),
                                                            landCover = rep(pop$cover, times = resIsla$samplingOccasions),
                                                            variable = rep(param, times = tr),
                                                            estimate = estimate, # As the growth rate is between periods, this is only from second on,
                                                            standardError = SE,
                                                            LCI = LCI,
                                                            UCI = UCI,
                                                            modelType = rep(resIsla$modelType, times = tr),
                                                            modelFormulation = rep(resIsla$bestModelName, times = tr),
                                                            AIC = rep(round(resIsla$bestModel@AIC, 2), times = tr),
                                                            K = rep(resIsla$bestModel@K, times = tr),
                                                            totalModelRuntime = rep(resIsla$totalModelRuntime, times = tr)
                                          )
                                          return(res)
                                          
                                        }))
        return(listParam)
      }))
      return(listIsland)
    }))
    write.csv(x = listSpecies, file = outputsTable)
  } else {
    message("Outputs table exists and rewriteTable = FALSE. Returning.")
    listSpecies <- data.table::fread(outputsTable, drop = "V1")
  }
  return(listSpecies)
}
