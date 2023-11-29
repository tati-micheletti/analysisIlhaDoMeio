calculatePopulationChange <- function(species, island = c("Ilha do Meio", 
                                                          "Ilha Rata"),
                                      envir = environment(),
                                      rewriteTable = FALSE){
  outputsTable <- file.path("./outputs/outputsTable.csv")
  if (any(!file.exists(outputsTable),
          rewriteTable)) {  
    spShort <- sapply(species, switch,
                       "Elaenia ridleyana" = "elaenia",
                       "Trachylepis atlantica" = "mabuia",
                       "Johngarthia lagostoma" = "crab",
                       "Sula dactylatra" = "maskedBooby")
  listSpecies <- rbindlist(lapply(spShort, function(sp){
    listIsland <- rbindlist(lapply(island, function(isla){
      islandShort <-  if (isla == "Ilha Rata") "Rata" else "Meio"
      resIsla <- get(paste0(sp, islandShort), envir = envir)
      listParam <- rbindlist(lapply(c("populationGrowthRate", 
                                      "detectionProbability",
                                      "immigration",
                                      "initialPopulationDensity"), function(param){
                                        message(paste0("Creating data for ", param, " for ", sp, " for ", isla))
                                        tr <- resIsla$samplingOccasions*resIsla$numberOfSites
                                        # POPULATION GROWTH
                                        if (param %in% c("populationGrowthRate", "initialPopulationDensity")){
                                          lng <- length(unique(resIsla$populationGrowthRate[,"Predicted"]))
                                          if (!any(lng == (resIsla$samplingOccasions-1), lng == 1)){
                                            message("Sites differ for population growth rate. Check.")
                                            browser()
                                          }
                                          popGrRate <- transpose(resIsla$populationGrowthRate[1:as.numeric(resIsla$samplingOccasions-1),])
                                          popGrRate <- cbind(rep(NA, NROW(popGrRate)), popGrRate)
                                          colnames(popGrRate) <- c(paste0("Cycle_", 1:resIsla$samplingOccasions))
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
                                          for (occ in 2:resIsla$samplingOccasions) {
                                            # Density
                                            newD <- newD*popGrRate["Predicted",paste0("Cycle_", occ)] # Make sure to check this is not constant through time...
                                            estimate <- c(estimate, newD)
                                            # SE
                                            newSE <- newSE*popGrRate["SE",paste0("Cycle_", occ)] # Make sure to check this is not constant through time...
                                            SE <- c(SE, newSE)
                                            # LCI
                                            newLCI <- newLCI*popGrRate["lower",paste0("Cycle_", occ)] # Make sure to check this is not constant through time...
                                            LCI <- c(LCI, newLCI)
                                            # UCI
                                            newUCI <- newUCI*popGrRate["upper",paste0("Cycle_", occ)] # Make sure to check this is not constant through time...
                                            UCI <- c(UCI, newUCI)
                                          }
                                          # Check that the length of newD matches the expected length
                                          if (length(estimate) != tr) 
                                            stop(paste0("Length of PopulationDensity doesn't match expected length: ",
                                                        tr,". Debug."))
                                          
                                        }
                                        # WE HAVE THE 40 SITES then the 5 CYCLES:
                                        # cycle 1: sites 1 to 40,
                                        # cycle 2: sites 1 to 40...
                                        res <- data.table(species = rep(resIsla$species, times = tr),
                                                          island = rep(isla, times = tr),
                                                          cycle = rep(1:resIsla$samplingOccasions, each = resIsla$numberOfSites),
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
    listSpecies <- data.table::fread(outputsTable)
  }
  return(listSpecies)
}
