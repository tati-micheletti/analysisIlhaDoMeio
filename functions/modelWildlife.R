############################### WILDLIFE MODELS ###############################
modelWildlife <- function(DT, island, species, rerunModels = FALSE, 
                              reRunK = FALSE){
  
  islandShort <-  if (island == "Ilha Rata") "Rata" else "Meio"
  spShort <-  if (species == "Elaenia ridleyana") "elaenia" else 
                if (species == "Trachylepis atlantica") "mabuia" else 
                  if (species == "Johngarthia lagostoma") "crab" else
                    if (species == "Sula dactylatra") "maskedBooby" else
                    stop(paste0("Species ", species, " is not in the dataset"))
  immg <- if (spShort %in% c("elaenia", "maskedBooby")) TRUE else FALSE
  paramTable <- loadParametersTable() 
  spParameters <- spParams(islandShort = islandShort, spShort = spShort)
  
  message(paste0("Running models for ", spShort, 
                 " on ", islandShort, "..."))
  
  ###### PREPARE DATA #####
  
  # First we subset the data to the interested place and the species
  if (spShort == "crab"){
    wildife <- DT[Species == species,]
    datesForMoon <- sort(unique(wildife$Date))
    # These were manually added: https://www.webexhibits.org/calendars/moon.html
    phasesOfMoon <- c("Waxing Gibbous Moon", "Full Moon", "Waning Crescent Moon",
                      "Waning Crescent Moon", "New Moon Crescent Moon", "New Moon Crescent Moon",
                      "Waxing Crescent Moon", "Waxing Gibbous Moon", "Waxing Gibbous Moon",
                      "Full Moon", "Waning Crescent Moon", 
                      "Waxing Crescent Moon")
    moon <- data.table(Date = as.Date(datesForMoon),
                       Moon = phasesOfMoon)
    wildifeMerged <- merge(wildife, moon, by = "Date")
    wildife <- wildifeMerged[Island == island,]
  } else {
    wildife <- DT[Island == island & Species == species,]
  }
  # If seabird species, use the whole island for observation as the whole island 
  # is walked for survey. Area here in meters, except for seabirds, which is m2. 
  area <- if (spShort == "maskedBooby") 180000 else 
                               unique(wildife$Radius)
  # Drop what we don't need, or don't want: 
  # Island, Year, Month, Day, Species, Radius, Original_Landscape, originDate
  toRemove <- c("Island", "Year", "Month", "Day", "Species", 
                "Radius", "Original_Landscape", "originDate") 
  wildife <- wildife[, (toRemove) := NULL]

  # Now we make the needed objects:
  # y --> Matrix with observations where rows represent each Site (unmarked calls this 'M'), and columns the 
  # (increasing) sampling occasion (i.e., the time series of counts, representing 
  # each visit), which unmarked calls 'T'. The value is the Counts.
  finalCounts <- dcast(wildife, Site ~ Date, value.var = "Counts")
  y <- as.matrix(finalCounts[,-1]) # Convert to matrix
  
  # obsCov --> Observation covariates are the ones that vary per sampling occasion.
  # This object is a named list of each covariate, which are data frames of rows 
  # representing each Site and the columns each sampling occasion. 
  # Examples are: observerID, JulianDate, TSE
  JulianDate <- as.matrix(dcast(wildife, Site ~ Date, value.var = "JulianDate")[,-1])
  TSE <- as.matrix(dcast(wildife, Site ~ Date, value.var = "timeSinceStartEradication")[,-1])
  
  # Fill in the matrix of JulianDate and TSE: this affects the models later on
  JulianDate <- replaceNAwithColMean(JulianDate)
  TSE <- replaceNAwithColMean(TSE)

  if (spShort %in% c("elaenia", "maskedBooby")) {
    obsCovs <- list(JulianDate = JulianDate,
                    TSE = TSE)
  } else {
    if (spShort == "crab"){
      moonPhase <- replaceNAwithColMean(as.matrix(dcast(wildife, Site ~ Date, value.var = "Moon")[,-1]))
      obsCovs <- list(JulianDate = JulianDate,
                      TSE = TSE,
                      moonPhase = moonPhase)
    } else { # Mabuia
        observerID <- as.matrix(dcast(wildife, Site ~ Date, value.var = "Observer")[,-1])
        obsCovs <- list(JulianDate = JulianDate,
                        TSE = TSE,
                        observerID = observerID)
    }
  }

  # siteCovs --> Site covariates are fixed through time, varying only among each 
  # other. This object is matrix where rows representing each Site and the columns 
  # each covariate. 
  # Examples are: Cover (landscape type)
  coverRaw <- dcast(wildife, Site ~ Date, value.var = "Landscape")
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
    numPrimary <- numPrimary/2
    # Since the closed periods are really close, we can use one variable for both days 
    # (yearlySiteCovs). To do this, we need to get the mean of every n columns (n = 
    # secondary periods, or 2 in our case). I will do this using a loop
    yearlySiteCovs <- list(JulianDate = as.matrix(joinCols(y = y, DT = JulianDate, 
                                                           numPrimary = numPrimary)),
                           TSE = as.matrix(joinCols(y = y, DT = TSE, 
                                                    numPrimary = numPrimary)))
  } else {
    if (spShort == "crab"){
    yearlySiteCovs <- list(moonPhase = moonPhase,
                           JulianDate = JulianDate,
                           TSE = TSE)
    } else {
      yearlySiteCovs <- list(JulianDate = JulianDate,
                             TSE = TSE)
    }
  }

  # Create the unmarkedFrame
  wildifeDF <- unmarkedFramePCO(y = y,
                               siteCovs = cover,
                               obsCovs = obsCovs,
                               numPrimary = numPrimary,
                               yearlySiteCovs = yearlySiteCovs
  )
  # Take a look
  wildifeDF
  summary(wildifeDF)
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
  
  wildife_nullModels_file <- paste0("./outputs/", spShort,"_", islandShort,"_nullModels.qs") 
  
  if (any(rerunModels, !file.exists(wildife_nullModels_file))){
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
                       data = wildifeDF,
                       mixture = nullModType, # Negative binomial
                       dynamics = "trend", # We want the population trend through time
                       immigration = immg,
                       K = K1)
      return(nM)
    })
    names(nullModels) <- nullModsNames
    qs::qsave(x = nullModels, file = wildife_nullModels_file)
  } else {
    message(paste0("Null models for ", spShort," on ", island,
                   " are available. Loading results..."))
    nullModels <- qs::qread(file = wildife_nullModels_file)
  }
  # Now we check which formulation to use. Runs all three formulations: ZIP, P, and NB
  # ZIP = Zero Inflated; NB = Negative Binomial; P = Poisson
  # ZIP is significant, Dispersion not: "ZIP"
  # ZIP is not significant, Dispersion is: "NB"
  # ZIP is not significant, Dispersion is not significant: "P"
  # ZIP is significant, Dispersion is significant: best AIC
  pZI2 <- summary(nullModels[["ZIP"]])[["psi"]][["P(>|z|)"]] < 0.05
  pNB <- summary(nullModels[["NB"]])[["alpha"]][["P(>|z|)"]] < 0.05
  
  bestAIC <- if (nullModels[["ZIP"]]@AIC < nullModels[["NB"]]@AIC) "ZIP" else "NB"
  
  modType <- if (all(any(pZI1, pZI2), !pNB)) "ZIP" else # ZIP = TRUE, NB = FALSE
    if (all(!any(pZI1, pZI2), pNB)) "NB" else # ZIP = FALSE, NB = TRUE
      if (all(!any(pZI1, pZI2), !pNB)) "P" else # ZIP = FALSE, NB = FALSE
        if (all(any(pZI1, pZI2), !pNB)) bestAIC else # ZIP = TRUE, NB = TRUE
          stop(paste0("NB = ", pNB, "; ZIP1 = ", pZI1, 
                      "; ZIP2 = ", pZI2, 
                      ". Something seems wrong. Debug."))
  message(paste0("Type of model for ", spShort, " for ", island, 
                 " was determined as ", modType))
  ######
  
  ###### MODELS #####
  wildife_models <- paste0("./outputs/",spShort,"_",islandShort,"_models.qs")
  wildife_models_t <- paste0("./outputs/",spShort,"_",islandShort,"_models_time.qs")
  
  mods <- data.table(expand.grid(lambda = spParameters$lambda, 
                                 gamma = spParameters$gamma, 
                                 p = spParameters$p, 
                                 iota = spParameters$iota))
  mods[, models := paste0(lambda, gamma, p, iota)]
  
  if (any(rerunModels, !file.exists(wildife_models))){
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
      indMod <- paste0("./outputs/individualModels/", spShort, "_", islandShort, "_", modName,".qs")
      if (!file.exists(indMod)){
        message(paste0("The model ", modName, " for ", spShort,
                       " for ", island, " doesn't exist. Running."))
        currMod <- pcountOpen(lambdaformula = paramTable[variable == as.character(mods[Row, lambda]), value][[1]],  # Initial abundance
                              gammaformula = paramTable[variable == as.character(mods[Row, gamma]), value][[1]],  # Population growth rate
                              omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                              pformula = paramTable[variable == as.character(mods[Row, p]), value][[1]], # Probability of observation
                              iotaformula = paramTable[variable == as.character(mods[Row, iota]), value][[1]], # Immigration
                              data = wildifeDF,
                              mixture = modType,
                              dynamics = "trend", # We want the population trend through time
                              immigration = immg,
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
    
    wildife_Models <- fitList(fits = allMods)
    
    wildife_models_time <- numeric(length(wildife_Models@fits))
    for (m in 0:(length(wildife_models_time)-1)){
      otm <- get(paste0("t", m+1))-get(paste0("t", m))
      tm <- as.numeric(otm)
      wildife_models_time[m+1] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    
    qs::qsave(x = wildife_Models, file = wildife_models)
    qs::qsave(x = wildife_models_time, file = wildife_models_t)
    
  } else {
    message(paste0("Models for ", spShort, " on ", island,
                   " are available. Loading results..."))
    wildife_Models <- qs::qread(file = wildife_models)
    wildife_models_time <- qs::qread(file = wildife_models_t)
  }
  
  message(paste0("It takes ", round(sum(wildife_models_time), 0), 
                 " minutes to run all models."))
  
  modelSelected <- modSel(wildife_Models)
  wildife_best_model_name <- modelSelected@Full[1,"model"]
  wildife_best_model <- wildife_Models@fits[[wildife_best_model_name]]
  
  ######
  
  ###### MODELS K #####
  
  # We need to check if our K is enough. This is done my running the best model 
  # with various values of K until AIC remains unchanged
  
  wildife_best <- paste0("./outputs/",spShort,"_",islandShort,"_bestModel.qs")
  wildife_models_K_t <- paste0("./outputs/",spShort,"_",islandShort,"_models_K_time.qs")
  # I need to chose K iteratively. The first K needs to be the one that is 
  # max observation + 1. Then, we can progress
  if (any(reRunK,
          rerunModels, 
          !file.exists(wildife_best))){
    if (any(reRunK, rerunModels))
      message(paste0("The argument rerunModels = TRUE or rerunK = TRUE. ", 
                     "Running the models for ", spShort," for ", island," for assessing K. It might take some time...")) else   
                       message(paste0("Models for assessing K for ", spShort, " on ",
                                      island, " haven't yet been run yet.",
                                      " Running the models will take some time..."))

    bestLambda <- as.character(mods[models == wildife_best_model_name, lambda])
    bestGamma <- as.character(mods[models == wildife_best_model_name, gamma])
    bestP <- as.character(mods[models == wildife_best_model_name, p])
    bestIota <- if (spShort == "elaenia") as.character(mods[models == wildife_best_model_name, iota]) else "iota(.)"
    bestModIndex <- which(names(wildife_Models@fits) == wildife_best_model_name)
    t0 <- wildife_models_time[bestModIndex] # Time it took the best model for first K
    t1 <- Sys.time()
      en <- environment()
      valuesK <- Kvals <- K1
      increase <- 200
      tolerance <- 0.1
      previousAIC <- -1
      currentAIC <- wildife_best_model@AIC
      i <- 2
      wildife_bestModel <- list()
      wildife_bestModel[[paste0("K_", Kvals)]] <- wildife_best_model
      while (!AICmatches(previousAIC, currentAIC, tolerance)){
        message(paste0("Current K is ", Kvals, ". Previous AIC was ", round(previousAIC, digits = 2),
                       ", while current AIC is ", round(currentAIC, digits = 2), ". With a tolerance ",
                       "of ", tolerance, " K is not yet stable. Trying K at ",
                       Kvals + increase))
        Kvals <- Kvals + increase
        valuesK <- c(valuesK, Kvals)
        previousAIC <- currentAIC
        tic(paste0("Model with K ", Kvals, " finished: "))
        indMod <- paste0("./outputs/individualModels/", spShort, "_", 
                         islandShort, "_K", Kvals,".qs")
        if (!file.exists(indMod)){
          wildife_bestModel[[paste0("K_", Kvals)]] <- pcountOpen(lambdaformula = paramTable[variable == bestLambda, value][[1]],  # Initial abundance
                           gammaformula = paramTable[variable == bestGamma, value][[1]],  # Formula for population growth rate
                           omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                           pformula = paramTable[variable == bestP, value][[1]],
                           iotaformula = paramTable[variable == bestIota, value][[1]],
                           data = wildifeDF,
                           mixture = modType,
                           dynamics = "trend", # We want the population trend through time
                           immigration = immg,
                           K = Kvals)
        toc()
        qs::qsave(x =  wildife_bestModel[[paste0("K_", Kvals)]], file = indMod)
        } else {
          message(paste0("The model with K ", Kvals, " for ", spShort,
                         " for ", island," already exists. Returning."))
          wildife_bestModel[[paste0("K_", Kvals)]] <- qs::qread(indMod)
        }
        # After making and saving the model, update the AIC
        currentAIC <- wildife_bestModel[[paste0("K_", Kvals)]]@AIC
        assign(x = paste0("t", i), value = Sys.time(), envir = en)
        i <- i + 1
      }
      message(paste0("Model stabilized at K ", Kvals, ". Previous AIC was ", 
                     round(previousAIC, digits = 2), ", while current AIC is ", 
                     round(currentAIC, digits = 2), "."))
    qs::qsave(x = wildife_bestModel, file = wildife_best)
    wildife_models_K_time <- numeric(length(wildife_bestModel))
    for (m in 2:length(wildife_bestModel)){
      otm <- get(paste0("t", m))-get(paste0("t", m-1))
      tm <- as.numeric(otm)
      wildife_models_K_time[m] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    wildife_models_K_time[1] <- t0
    qs::qsave(x = wildife_models_K_time, file = wildife_models_K_t)
  } else {
    message(paste0("Models for assessing K for ", spShort," on ", island,
                   " are available. Loading results..."))
    wildife_bestModel <- qs::qread(file = wildife_best)
    wildife_models_K_time <- qs::qread(file = wildife_models_K_t)
  }
  wildife_AICs <- sapply(wildife_bestModel, function(m) m@AIC)
  
  # Plot AIC to see if the K is ok
  p <- plot(x = valuesK, y = wildife_AICs, xlab = "K values", ylab = "AIC", 
            title = paste0(spShort, " models on ", island))

  # We then override the best model with the model with highest K --> most stable parameters
  wildife_best_model <- wildife_bestModel[[length(wildife_bestModel)]]
  ######
  
  ###### PARAMETER ESTIMATION ######
  # Back-transformation of parameters using predict()
  initialPopulation <- predict(obj = wildife_best_model, type = "lambda")
  populationGrowthRate <- predict(wildife_best_model, type = "gamma")
  detection <- predict(wildife_best_model, type = "det")
  immigration <- if (spShort == "elaenia") predict(wildife_best_model, type = "iota") else NA
  
  pop <- unique(data.table(initialPop = predict(wildife_best_model, type = "lambda")[,1], 
                           cover = wildife_best_model@data@siteCovs[, "cover"]))
  
  # The observation radius for wildife is 3m, so the area is 28.27m2 or 0,002827ha
  if (spShort %in% c("elaenia", "maskedBooby")){
    pop[, areaM := area] # It is already correct for boobies (i.e., whole island)
  } else {
    pop[, areaM := (1*pi*area^2)] # height*pi*observation radius^2 = observed area --> 
    # Observations for mabuia were done in 3D, trees, rocks, bushes, etc.
  }
  pop[, densityM := initialPop/areaM]
  ######
  
  return(list(IndividualModelRunTime = wildife_models_time,
              totalModelRunTime = sum(wildife_models_time),
              nullModels = nullModels,
              allModels = wildife_Models,
              bestModelName = wildife_best_model_name,
              bestModel = wildife_best_model,
              allKmodels = wildife_bestModel,
              Kplot = p,
              initialPopulation = initialPopulation,
              detectionProbability = detection,
              populationGrowthRate = populationGrowthRate,
              immigration = immigration,
              initialPopulationDensity = pop))
}
#################################
