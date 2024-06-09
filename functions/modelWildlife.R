############################### WILDLIFE MODELS ###############################
modelWildlife <- function(DT, island, species, rerunModels = FALSE, 
                              reRunK = FALSE, plotting = FALSE,
                          weatherData, useImmigration,
                          useRobustDesign){
  
  ###### GENERAL SETUP #####
  
  islandShort <-  if (island == "Ilha Rata") "Rata" else "Meio"
  spShort <-  if (species == "Elaenia ridleyana") "elaenia" else 
                if (species == "Trachylepis atlantica") "mabuia" else 
                  if (species == "Johngarthia lagostoma") "crab" else
                    if (species == "Sula dactylatra") "maskedBooby" else
                    stop(paste0("Species ", species, " is not in the dataset"))

  imm <- if (useImmigration) spShort else NULL
  paramTable <- loadParametersTable() 
  spParameters <- spParams(islandShort = islandShort, spShort = spShort, 
                           addImmigrationFor = imm)
  
  message(paste0("Running models for ", spShort, 
                 " on ", islandShort, "..."))
  
  toleranceSp <- getTolerance(sp = spShort,
                              isld = islandShort)
 
  addOns <- if (all(useImmigration, useRobustDesign)) "Immig_RbstDsgn" else 
    if (all(useRobustDesign, !useImmigration)) "RbstDsgn" else 
      if (all(!useRobustDesign, useImmigration)) "Immig" else NULL
  individualModels <- paste("individualModels", addOns, sep = "_")
  
  ###### PREPARE DATA #####
  
  # First we subset the data to the interested place and the species
  # [UPDATE]: Now I am using real data from the metheorological station of the island. 
  # Not necessary anymore to hardcode this below 
  # if (spShort == "crab"){
  #   wildlife <- DT[Species == species,]
  #   datesForMoon <- sort(unique(wildlife$Date))
  #   # These were manually added: https://www.webexhibits.org/calendars/moon.html
  #   phasesOfMoon <- c("Waxing Gibbous Moon", "Full Moon", "Waning Crescent Moon",
  #                     "Waning Crescent Moon", "New Moon Crescent Moon", "New Moon Crescent Moon",
  #                     "Waxing Crescent Moon", "Waxing Gibbous Moon", "Waxing Gibbous Moon",
  #                     "Full Moon", "Waning Crescent Moon", 
  #                     "Waxing Crescent Moon")
  #   moon <- data.table(Date = as.Date(datesForMoon),
  #                      Moon = phasesOfMoon)
  #   wildlifeMerged <- merge(wildlife, moon, by = "Date")
  #   wildlife <- wildlifeMerged[Island == island,]
  # } else {
    wildlife <- DT[Island == island & Species == species,]
    # }
  
  # Collate DT with weatherData
  ### Create totalRainfallThreeMonths from precip
  # 1. Create in the wildlife the column "startDate", which is 3 months from the "Date"
  wildlife[, startDate := Date-91]
  allDates <- weatherData$Date
  # 2. For each row in wildlife: 
  for (row in 1:NROW(wildlife)){
  #  2.1. Subset on weatherData the values for precip in which the "weatherData$Date" 
    #       is between "wildlife$startDate" and "wildlife$Date"
    currentDate <- wildlife[row, Date]
    startDate <- wildlife[row, startDate]
    endDate <- wildlife[row, Date]-1
    selectedDates <- subset(allDates, allDates >= startDate & allDates <= endDate)
    prepVals <- weatherData[Date %in% selectedDates, precip]
  #  2.2. Sum all values found in 2.1.
  #  2.3. Fill in the new column totalRainfallThreeMonths in wildlife with it 
    wildlife[row, totalRainfallThreeMonths := sum(prepVals)]
  ### Fill in moonphase
  # 2.1. While at each row to do the previous, match "wildlife$Date" and "weatherData$Date" 
  #      getting `moonphase` and fill the new `moonphase` variable in wildlife with the 
  #      value from weatherData  
    wildlife[row, moonphase := weatherData[Date == currentDate, moonphase]]
  }  
# If seabird species, use the whole island for observation as the whole island 
  # is walked for survey. Area here in meters, except for seabirds, which is m2. 
  area <- if (spShort == "maskedBooby") 180000 else 
                               unique(wildlife$Radius)
  # Drop what we don't need, or don't want: 
  # Island, Year, Month, Day, Species, Radius, Original_Landscape, originDate
  toRemove <- c("Island", "Year", "Month", "Day", "Species", "startDate",
                "Radius", "Original_Landscape", "originDate", "JulianDate") 
  wildlife <- wildlife[, (toRemove) := NULL]
  
  # Now we make the needed objects:
  # y --> Matrix with observations where rows represent each Site (unmarked calls this 'M'), and columns the 
  # (increasing) sampling occasion (i.e., the time series of counts, representing 
  # each visit), which unmarked calls 'T'. The value is the Counts.
  finalCounts <- dcast(wildlife, Site ~ Date, value.var = "Counts")
  y <- as.matrix(finalCounts[,-1]) # Convert to matrix
  
  # obsCov --> Observation covariates are the ones that vary per sampling occasion.
  # This object is a named list of each covariate, which are data frames of rows 
  # representing each Site and the columns each sampling occasion. 
  # Examples are: observerID, totalRainfallThreeMonths, TSE
  totalRainfallThreeMonths <- as.matrix(dcast(wildlife, Site ~ Date, value.var = "totalRainfallThreeMonths")[,-1])
  TSE <- as.matrix(dcast(wildlife, Site ~ Date, value.var = "timeSinceStartEradication")[,-1])
  
  # Fill in the matrix of totalRainfallThreeMonths and TSE: this affects the models later on
  totalRainfallThreeMonths <- replaceNAwithColMean(totalRainfallThreeMonths)
  TSE <- replaceNAwithColMean(TSE)

  if (spShort %in% c("elaenia", "maskedBooby")) {
    obsCovs <- list(totalRainfallThreeMonths = totalRainfallThreeMonths,
                    TSE = TSE)
  } else {
    if (spShort == "crab"){
      moonPhase <- replaceNAwithColMean(as.matrix(dcast(wildlife, Site ~ Date, value.var = "moonphase")[,-1]))
      obsCovs <- list(totalRainfallThreeMonths = totalRainfallThreeMonths,
                      TSE = TSE,
                      moonPhase = moonPhase)
    } else { # Mabuia
        observerID <- as.matrix(dcast(wildlife, Site ~ Date, value.var = "Observer")[,-1])
        obsCovs <- list(totalRainfallThreeMonths = totalRainfallThreeMonths,
                        TSE = TSE,
                        observerID = observerID)
    }
  }

  # siteCovs --> Site covariates are fixed through time, varying only among each 
  # other. This object is matrix where rows representing each Site and the columns 
  # each covariate. 
  # Examples are: Cover (landscape type)
  coverRaw <- dcast(wildlife, Site ~ Date, value.var = "Landscape")
  # For each row I need the unique value of cover without NA:
  cover <- data.frame(cover = unlist(lapply(1:NROW(coverRaw), 
                                            function(r){
                                              unique(na.omit(as.character(coverRaw[r,2:NCOL(coverRaw)])))
                                            })))
  # numPrimary --> Number of primary periods
  numPrimary <- NCOL(y)
  if (spShort == "elaenia"){
    # In the case of Elenia, we have 2 counts (i.e., secondary 
    # periods) per primary occasion
    # Since the closed periods are really close, we can use one variable for both days 
    # (yearlySiteCovs). To do this, we need to get the mean of every n columns (n =
    # secondary periods, or 2 in our case). I will do this using a loop
    
    # ROBUST DESIGN
    if (useRobustDesign){
      numPrimary <- numPrimary/2
      yearlySiteCovs <- list(totalRainfallThreeMonths = as.matrix(joinCols(y = y, DT = totalRainfallThreeMonths,
                                                             numPrimary = numPrimary)),
                             TSE = as.matrix(joinCols(y = y, DT = TSE,
                                                      numPrimary = numPrimary)))
    } else {
      yearlySiteCovs <- list(totalRainfallThreeMonths = totalRainfallThreeMonths,
                             TSE = TSE)
    }
  } else {
    if (spShort == "crab"){
    yearlySiteCovs <- list(moonPhase = moonPhase,
                           totalRainfallThreeMonths = totalRainfallThreeMonths,
                           TSE = TSE)
    } else {
      yearlySiteCovs <- list(totalRainfallThreeMonths = totalRainfallThreeMonths,
                             TSE = TSE)
    }
  }

  # Create the unmarkedFrame
  wildlifeDF <- unmarkedFramePCO(y = y,
                               siteCovs = cover,
                               obsCovs = obsCovs,
                               numPrimary = numPrimary,
                               yearlySiteCovs = yearlySiteCovs
  )
  # Take a look
  wildlifeDF
  summary(wildlifeDF)
  ######
  
  ###### NULL MODELS #####
  # Here we test for Zero Inflation on the data
  pZI1 <- test_ZI(dts = na.omit(as.numeric(y)))
  
  # We also need to chose K. K needs to be bigger than the max number observed, 
  # and it needs to be chosed as such as the increase in K doens't change the 
  # AIC's anymore. Here we go from max number + 1 or 50, whatever is bigger
  K1 <- max(max(y, na.rm = TRUE)+1, 50)
  message(paste0("K was defined as ", K1, " for ", spShort,
                 " for ", island))
  # We run then a negative binomial (NB) and a zero inflated poisson (ZIP)
  # We don't run Poisson because the data us zero inflated
  
  wildlife_nullModels_file <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                                                create = TRUE), 
                                        paste0(spShort,"_", islandShort,"_nullModels.qs"))
  
  if (any(rerunModels, !file.exists(wildlife_nullModels_file))){
    if (rerunModels)
      message(paste0("The argument rerunModels = TRUE.", 
                     "Running the null models for ", spShort," on ", island,
                     " should not take too much time...")) else   
                       message(paste0("Null models for ", spShort," on ", island,
                                      " haven't yet been run yet.",
                                      " Running these models should not take too much time..."))
    t0 <- Sys.time()
    nullModsNames <- c("NB", "ZIP", "P")
    nullModels <- lapply(nullModsNames, function(nullModType){
      nM <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~1,  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~1,
                       data = wildlifeDF,
                       mixture = nullModType, # Negative binomial
                       dynamics = "trend", # We want the population trend through time
                       immigration = useImmigration,
                       K = K1)
      return(nM)
    })
    names(nullModels) <- nullModsNames
    qs::qsave(x = nullModels, file = wildlife_nullModels_file)
  } else {
    message(paste0("Null models for ", spShort," on ", island,
                   " are available. Loading results..."))
    nullModels <- qs::qread(file = wildlife_nullModels_file)
  }
  # Now we check which formulation to use. Runs all three formulations: ZIP, P, and NB
  # ZIP = Zero Inflated; NB = Negative Binomial; P = Poisson
  # ZIP is significant, Dispersion not: "ZIP"
  # ZIP is not significant, Dispersion is: "NB"
  # ZIP is not significant, Dispersion is not significant: "P"
  # ZIP is significant, Dispersion is significant: best AIC
  pZI2 <- summary(nullModels[["ZIP"]])[["psi"]][["P(>|z|)"]] < 0.05
  pNB <- summary(nullModels[["NB"]])[["alpha"]][["P(>|z|)"]] < 0.05
  if (is.na(pZI2)) pZI2 <- FALSE
  if (is.na(pNB)) pNB <- FALSE
  bestAIC <- if (nullModels[["ZIP"]]@AIC < nullModels[["NB"]]@AIC) "ZIP" else "NB"

  modType <- if (all(!any(pZI1, pZI2), !pNB)) "P" else # ZIP = FALSE, NB = FALSE
    if (all(!any(pZI1, pZI2), pNB)) "NB" else # ZIP = FALSE, NB = TRUE
      if (all(any(pZI1, pZI2), !pNB)) "ZIP" else # ZIP = TRUE, NB = FALSE
        if (all(any(pZI1, pZI2), pNB)) bestAIC else # ZIP = TRUE, NB = TRUE
          stop(paste0("NB = ", pNB, "; ZIP1 = ", pZI1, 
                      "; ZIP2 = ", pZI2, 
                      ". Something seems wrong. Debug."))

  nullModChosen <- nullModels[[modType]]
  hasNoEstimates <- coef(nullModChosen)["gamTrend(Int)"] == 0
  if (hasNoEstimates){
    message(paste0("Model type ", modType, " did not converge. Using ", 
                   bestAIC, " instead"))
    modType <- bestAIC
  }
    message(paste0("Type of model for ", spShort, " for ", island, 
                   " was determined as ", modType))
  ######
  
  ###### MODELS #####
  wildlife_models <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                                       create = TRUE), paste0(spShort,"_",islandShort,"_models.qs"))
  wildlife_models_t <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                                         create = TRUE), paste0(spShort,"_",islandShort,"_models_time.qs"))
  
  mods <- data.table(expand.grid(lambda = spParameters$lambda, 
                                 gamma = spParameters$gamma, 
                                 p = spParameters$p, 
                                 iota = spParameters$iota))
  mods[, models := paste0(lambda, gamma, p, iota)]
  
  if (any(rerunModels, !file.exists(wildlife_models))){
    if (rerunModels)
      message(paste0("The argument rerunModels = TRUE.", 
                     "Running the models will take some time...")) else   
                       message(paste0("Models for ", spShort," on ", island,
                                      " haven't yet been run yet.",
                                      " Running the models will take some time..."))
    t0 <- Sys.time()
    en <- environment()
    allMods <- lapply(1:NROW(mods), function(Row){
      modName <- as.character(mods[Row, models])
      tic(paste0("Model ", modName, " finished: "))
      indMod <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                        create = TRUE), paste0(spShort, "_", islandShort, "_", modName,".qs"))
      if (any(!file.exists(indMod),
              rerunModels)){
        message(paste0("The model ", modName, " for ", spShort,
                       " for ", island, " doesn't exist. Running."))
        currMod <- pcountOpen(lambdaformula = paramTable[variable == as.character(mods[Row, lambda]), value][[1]],  # Initial abundance
                              gammaformula = paramTable[variable == as.character(mods[Row, gamma]), value][[1]],  # Population growth rate
                              omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                              pformula = paramTable[variable == as.character(mods[Row, p]), value][[1]], # Probability of observation
                              iotaformula = paramTable[variable == as.character(mods[Row, iota]), value][[1]], # Immigration
                              data = wildlifeDF,
                              mixture = modType,
                              dynamics = "trend", # We want the population trend through time
                              immigration = useImmigration,
                              K = K1)
        toc()
        qs::qsave(x = currMod, file = indMod)
      } else {
        message(paste0("The model ", modName, " already exists. Returning."))
        currMod <- qs::qread(indMod)
      }
      assign(x = paste0("t", Row), value = Sys.time(), envir = en)
      return(currMod)
    })
    names(allMods) <- mods[["models"]]
    
    wildlife_Models <- fitList(fits = allMods)
    
    wildlife_models_time <- numeric(length(wildlife_Models@fits))
    for (m in 0:(length(wildlife_models_time)-1)){
      otm <- get(paste0("t", m+1))-get(paste0("t", m))
      tm <- as.numeric(otm)
      wildlife_models_time[m+1] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    
    qs::qsave(x = wildlife_Models, file = wildlife_models)
    qs::qsave(x = wildlife_models_time, file = wildlife_models_t)
    
  } else {
    message(paste0("Models for ", spShort, " on ", island,
                   " are available. Loading results..."))
    wildlife_Models <- qs::qread(file = wildlife_models)
    wildlife_models_time <- qs::qread(file = wildlife_models_t)
  }
  
  message(paste0("It takes ", round(sum(wildlife_models_time), 0), 
                 " minutes to run all models."))
  
  modelSelected <- modSel(wildlife_Models)
  wildlife_best_model_name <- modelSelected@Full[1,"model"]
  wildlife_best_model <- wildlife_Models@fits[[wildlife_best_model_name]]
  
  ######
  
  ###### MODELS K #####
  
  # We need to check if our K is enough. This is done my running the best model 
  # with various values of K until AIC remains unchanged
  wildlife_best <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                                     create = TRUE), paste0(spShort,"_",islandShort,"_bestModel.qs"))
  wildlife_models_K_t <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), 
                                                           create = TRUE), paste0(spShort,"_",islandShort,"_models_K_time.qs"))
  # I need to chose K iteratively. The first K needs to be the one that is 
  # max observation + 1. Then, we can progress
  if (any(reRunK,
          rerunModels, 
          !file.exists(wildlife_best))){
    if (any(reRunK, rerunModels))
      message(paste0("The argument rerunModels = TRUE or rerunK = TRUE. ", 
                     "Running the models for ", spShort," for ", island," for assessing K. It might take some time...")) else   
                       message(paste0("Models for assessing K for ", spShort, " on ",
                                      island, " haven't yet been run yet.",
                                      " Running the models will take some time..."))

    bestLambda <- as.character(mods[models == wildlife_best_model_name, lambda])
    bestGamma <- as.character(mods[models == wildlife_best_model_name, gamma])
    bestP <- as.character(mods[models == wildlife_best_model_name, p])
    bestIota <- if (useImmigration) as.character(mods[models == wildlife_best_model_name, iota]) else "iota(.)"
    bestModIndex <- which(names(wildlife_Models@fits) == wildlife_best_model_name)
    t0 <- wildlife_models_time[bestModIndex] # Time it took the best model for first K
    t1 <- Sys.time()
      en <- environment()
      valuesK <- Kvals <- K1
      increase <- 200
      cmbSpIs <- paste0(spShort, "_", islandShort)
      tolerance <- toleranceSp
      previousAIC <- -1
      currentAIC <- wildlife_best_model@AIC
      i <- 2
      wildlife_bestModel <- list()
      wildlife_bestModel[[paste0("K_", Kvals)]] <- wildlife_best_model
      while (!AICmatches(previousAIC, currentAIC, tolerance)){
        message(paste0("Current K is ", Kvals, ". Previous AIC was ", round(previousAIC, digits = 2),
                       ", while current AIC is ", round(currentAIC, digits = 2), ". With a tolerance ",
                       "of ", tolerance, " K is not yet stable. Trying K at ",
                       Kvals + increase))
        Kvals <- Kvals + increase
        valuesK <- c(valuesK, Kvals)
        previousAIC <- currentAIC
        tic(paste0("Model with K ", Kvals, " finished: "))
        
        indMod <- file.path(reproducible::checkPath(file.path("./outputs", individualModels), create = TRUE), 
                            paste0(spShort, "_", islandShort, "_K", Kvals,".qs"))
        if (any(!file.exists(indMod),
                reRunK, rerunModels)){
          wildlife_bestModel[[paste0("K_", Kvals)]] <- pcountOpen(lambdaformula = paramTable[variable == bestLambda, value][[1]],  # Initial abundance
                           gammaformula = paramTable[variable == bestGamma, value][[1]],  # Formula for population growth rate
                           omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                           pformula = paramTable[variable == bestP, value][[1]],
                           iotaformula = paramTable[variable == bestIota, value][[1]],
                           data = wildlifeDF,
                           mixture = modType,
                           dynamics = "trend", # We want the population trend through time
                           immigration = useImmigration,
                           K = Kvals)
        toc()
        qs::qsave(x =  wildlife_bestModel[[paste0("K_", Kvals)]], file = indMod)
        } else {
          message(paste0("The model with K ", Kvals, " for ", spShort,
                         " for ", island," already exists. Returning."))
          wildlife_bestModel[[paste0("K_", Kvals)]] <- qs::qread(indMod)
        }
        # After making and saving the model, update the AIC
        currentAIC <- wildlife_bestModel[[paste0("K_", Kvals)]]@AIC
        assign(x = paste0("t", i), value = Sys.time(), envir = en)
        i <- i + 1
      }
      message(paste0("Model stabilized at K ", Kvals, ". Previous AIC was ", 
                     round(previousAIC, digits = 2), ", while current AIC is ", 
                     round(currentAIC, digits = 2), "."))
    qs::qsave(x = wildlife_bestModel, file = wildlife_best)
    wildlife_models_K_time <- numeric(length(wildlife_bestModel))
    for (m in 2:length(wildlife_bestModel)){
      otm <- get(paste0("t", m))-get(paste0("t", m-1))
      tm <- as.numeric(otm)
      wildlife_models_K_time[m] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    wildlife_models_K_time[1] <- t0
    qs::qsave(x = wildlife_models_K_time, file = wildlife_models_K_t)
  } else {
    message(paste0("Models for assessing K for ", spShort," on ", island,
                   " are available. Loading results..."))
    wildlife_bestModel <- qs::qread(file = wildlife_best)
    wildlife_models_K_time <- qs::qread(file = wildlife_models_K_t)
  }
  wildlife_AICs <- sapply(wildlife_bestModel, function(m) m@AIC)
  
  # Plot AIC to see if the K is ok
  if (plotting){
    p <- plot(x = valuesK, y = wildlife_AICs, xlab = "K values", ylab = "AIC", 
              title = paste0(spShort, " models on ", island))
  } else p <- NA

  # We then override the best model with the model with highest K --> most stable parameters
  wildlife_best_model <- wildlife_bestModel[[length(wildlife_bestModel)]]
  ######
  
  ###### PARAMETER ESTIMATION ######
  # Back-transformation of parameters using predict()
  initialPopulation <- predict(obj = wildlife_best_model, type = "lambda")
  populationGrowthRate <- predict(wildlife_best_model, type = "gamma")
  detection <- predict(wildlife_best_model, type = "det")
  if (useImmigration){
    immigration <- predict(wildlife_best_model, type = "iota")
  } else immigration <- NA 
  
  pop <- unique(data.table(initialPop = predict(wildlife_best_model, type = "lambda")[,1], 
                           cover = wildlife_best_model@data@siteCovs[, "cover"]))
  
  # The observation radius for wildlife is 3m, so the area is 28.27m2 or 0,002827ha
  if (spShort %in% c("elaenia", "maskedBooby")){
    pop[, areaM := area] # It is already correct for boobies (i.e., whole island)
  } else {
    pop[, areaM := (1*pi*area^2)] # height*pi*observation radius^2 = observed area --> 
    # Observations for mabuia were done in 3D, trees, rocks, bushes, etc.
  }
  pop[, densityM := initialPop/areaM]
  ######
  
  return(list(species = species,
              island = island,
              nullModels = nullModels,
              allModels = wildlife_Models,
              bestModelName = wildlife_best_model_name,
              bestModel = wildlife_best_model,
              allKmodels = wildlife_bestModel,
              Kplot = p,
              initialPopulation = initialPopulation,
              detectionProbability = detection,
              populationGrowthRate = populationGrowthRate,
              immigration = immigration,
              initialPopulationDensity = pop,
              modelType = modType,
              numberOfSites = NROW(y),
              cover = wildlife_best_model@data@siteCovs[, "cover"],
              samplingOccasions = numPrimary,
              totalModelRuntime = sum(wildlife_models_K_time[2:length(wildlife_models_K_time)], 
                                          wildlife_models_time, na.rm = TRUE)))
}
#################################
