############################### WILDLIFE MODELS ###############################
modelTerrestrials <- function(DT, island, species, rerunModels = FALSE, reRunK = FALSE){
  
  islandShort <-  if (island == "Ilha Rata") "Rata" else "Meio"
  spShort <-  if (species == "Elaenia ridleyana") "elaenia" else 
                if (species == "Trachylepis atlantica") "mabuia" else 
                  if (species == "Johngarthia lagostoma") "crab" else 
                    stop(paste0("Species ", species, " is not in the dataset"))
  immg <- if (species == "Elaenia ridleyana") TRUE else FALSE
  paramTable <- loadParametersTable() 
  spParameters <- spParams(islandShort = islandShort, spShort = spShort)
  
  message(paste0("Running models for ", spShort, 
                 " on ", islandShort, "..."))
  
  ###### PREPARE DATA #####
  
  # First we subset the data to the interested place and the species
  if (spShort == "crab"){
    terrestrials <- DT[Species == species,]
    datesForMoon <- sort(unique(terrestrials$Date))
    # These were manually added: https://www.webexhibits.org/calendars/moon.html
    phasesOfMoon <- c("Waxing Gibbous Moon", "Full Moon", "Waning Crescent Moon",
                      "Waning Crescent Moon", "New Moon Crescent Moon", "New Moon Crescent Moon",
                      "Waxing Crescent Moon", "Waxing Gibbous Moon", "Waxing Gibbous Moon",
                      "Full Moon", "Waning Crescent Moon", "Waning Crescent Moon", 
                      "Waxing Crescent Moon")
    moon <- data.table(Date = as.Date(datesForMoon),
                       Moon = phasesOfMoon)
    terrestrialsMerged <- merge(terrestrials, moon, by = "Date")
    terrestrials <- terrestrialsMerged[Island == island,]
  } else {
    terrestrials <- DT[Island == island & Species == species,]
  }
  area <- unique(terrestrials$Radius)
  # Drop what we don't need, or don't want: 
  # Island, Year, Month, Day, Species, Radius, Original_Landscape, originDate
  toRemove <- c("Island", "Year", "Month", "Day", "Species", 
                "Radius", "Original_Landscape", "originDate",
                "pointID", "transectID") # For terrestrials in Rata, we need to use ponto_numero as site
  # because there were mistakes in the table that were corrected
  # by PM later.
  if (islandShort == "Rata"){
    toRemove <- c(toRemove, "site")
    terrestrials <- terrestrials[, (toRemove) := NULL]
    names(terrestrials)[names(terrestrials)=="ponto_numero"] <- "site"
    terrestrials[, site := as.character(site)]
  } else {
    toRemove <- c(toRemove, "ponto_numero")
    terrestrials <- terrestrials[, (toRemove) := NULL]
  }
  
  # Now we make the needed objects:
  # y --> Matrix with observations where rows represent each site (unmarked calls this 'M'), and columns the 
  # (increasing) sampling occasion (i.e., the time series of counts, representing 
  # each visit), which unmarked calls 'T'. The value is the Counts.
  finalCounts <- dcast(terrestrials, site ~ Date, value.var = "Counts")
  y <- as.matrix(finalCounts[,-1]) # Convert to matrix
  
  # obsCov --> Observation covariates are the ones that vary per sampling occasion.
  # This object is a named list of each covariate, which are data frames of rows 
  # representing each site and the columns each sampling occasion. 
  # Examples are: observerID, JulianDate, TSE
  moonPhase <- if (spShort == "crab") as.matrix(dcast(terrestrials, site ~ Date, value.var = "Moon")[,-1]) else NA
  observerID <- as.matrix(dcast(terrestrials, site ~ Date, value.var = "Observer")[,-1])
  JulianDate <- as.matrix(dcast(terrestrials, site ~ Date, value.var = "JulianDate")[,-1])
  TSE <- as.matrix(dcast(terrestrials, site ~ Date, value.var = "timeSinceStartEradication")[,-1])
  
  # Fill in the matrix of JulianDate and TSE: this affects the models later on
  JulianDate <- replaceNAwithColMean(JulianDate)
  TSE <- replaceNAwithColMean(TSE)
  if (spShort == "crab") moonPhase <- replaceNAwithColMean(moonPhase)
  
  if (spShort == "elaenia") {
    obsCovs <- list(JulianDate = JulianDate,
                    TSE = TSE)
  } else {
    if (spShort == "crab"){
      obsCovs <- list(JulianDate = JulianDate,
                      TSE = TSE,
                      moonPhase = moonPhase)
    } else {
      obsCovs <- list(JulianDate = JulianDate,
                      TSE = TSE,
                      observerID = observerID)
    }
  }

  # siteCovs --> Site covariates are fixed through time, varying only among each 
  # other. This object is matrix where rows representing each site and the columns 
  # each covariate. 
  # Examples are: Cover (landscape type)
  coverRaw <- dcast(terrestrials, site ~ Date, value.var = "Landscape")
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
    interval <- NCOL(y)/numPrimary # How many secondary periods we have
    firstColsToJoin <- seq(from = 1, to = NCOL(y), by = interval)
    # For TSE
    newTSE <- lapply(firstColsToJoin, function(i){
      colsToJoin <- colnames(TSE)[i:(i+1)]
      TSE <- data.table(TSE[, colsToJoin])
      TSE[, colnames(TSE)[i] := rowMeans(.SD), .SDcols = names(TSE)]
      TSE[,names(TSE)[NCOL(TSE)] := NULL]
      return(TSE)
    })
    newTSE <- setDT(unlist(newTSE, recursive = FALSE),check.names = TRUE)[]
    # For JulianDate
    newJD <- lapply(firstColsToJoin, function(i){
      colsToJoin <- colnames(JulianDate)[i:(i+1)]
      JD <- data.table(JulianDate[, colsToJoin])
      JD[, colnames(JulianDate)[i] := rowMeans(.SD), .SDcols = names(JD)]
      JD[,names(JD)[NCOL(JD)] := NULL]
      return(JD)
    })
    newJD <- setDT(unlist(newJD, recursive = FALSE),check.names = TRUE)[]
    
    yearlySiteCovs <- list(JulianDate = as.matrix(newJD),
                           TSE = as.matrix(newTSE))
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
  terrestrialsDF <- unmarkedFramePCO(y = y,
                               siteCovs = cover,
                               obsCovs = obsCovs,
                               numPrimary = numPrimary,
                               yearlySiteCovs = yearlySiteCovs
  )
  # Take a look
  terrestrialsDF
  summary(terrestrialsDF)
  
  ######
  
  ###### NULL MODELS #####
  # Here we test for Zero Inflation on the data
  pZI1 <- test_ZI(dts = na.omit(as.numeric(y)))
  
  # We run then a negative binomial (NB) and a zero inflated poisson (ZIP)
  # We don't run Poisson because the data us zero inflated
  
  terrestrials_nullModels_file <- paste0("./outputs/", spShort,"_", islandShort,"_nullModels.qs") 
  
  if (any(rerunModels, !file.exists(terrestrials_nullModels_file))){
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
                       data = terrestrialsDF,
                       mixture = nullModType, # Negative binomial
                       dynamics = "trend", # We want the population trend through time
                       immigration = immg,
                       K = 100)
      return(nM)
    })
    names(nullModels) <- nullModsNames
    qs::qsave(x = nullModels, file = terrestrials_nullModels_file)
  } else {
    message(paste0("Null models for terrestrials on ", island," are available. Loading results..."))
    nullModels <- qs::qread(file = terrestrials_nullModels_file)
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
  terrestrials_models <- paste0("./outputs/",spShort,"_",islandShort,"_models.qs")
  terrestrials_models_t <- paste0("./outputs/",spShort,"_",islandShort,"_models_time.qs")
  
  mods <- data.table(expand.grid(lambda = spParameters$lambda, 
                                 gamma = spParameters$gamma, 
                                 p = spParameters$p, 
                                 iota = spParameters$iota))
  mods[, models := paste0(lambda, gamma, p, iota)]
  
  if (any(rerunModels, !file.exists(terrestrials_models))){
    if (rerunModels)
      message(paste0("The argument rerunModels = TRUE.", 
                     "Running the models will take some time...")) else   
                       message(paste0("Models for terrestrials on ", island," haven't yet been run yet.",
                                      " Running the models will take some time..."))
    t0 <- Sys.time()
    en <- environment()
    allMods <- lapply(1:NROW(mods), function(Row){
      modName <- as.character(mods[Row, models])
      tic(paste0("Model ", modName, " finished: "))
      indMod <- paste0("./outputs/individualModels/", spShort, "_", islandShort, "_", modName,".qs")
      if (!file.exists(indMod)){
        message(paste0("The model ", modName, " doesn't exist. Running."))
        currMod <- pcountOpen(lambdaformula = paramTable[variable == as.character(mods[Row, lambda]), value][[1]],  # Initial abundance
                              gammaformula = paramTable[variable == as.character(mods[Row, gamma]), value][[1]],  # Population growth rate
                              omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                              pformula = paramTable[variable == as.character(mods[Row, p]), value][[1]], # Probability of observation
                              iotaformula = paramTable[variable == as.character(mods[Row, iota]), value][[1]], # Immigration
                              data = terrestrialsDF,
                              mixture = modType,
                              dynamics = "trend", # We want the population trend through time
                              immigration = immg,
                              K = 100)
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
    
    terrestrials_Models <- fitList(fits = allMods)
    
    terrestrials_models_time <- numeric(length(terrestrials_Models@fits))
    for (m in 0:(length(terrestrials_models_time)-1)){
      otm <- get(paste0("t", m+1))-get(paste0("t", m))
      tm <- as.numeric(otm)
      terrestrials_models_time[m+1] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    
    qs::qsave(x = terrestrials_Models, file = terrestrials_models)
    qs::qsave(x = terrestrials_models_time, file = terrestrials_models_t)
    
  } else {
    message(paste0("Models for ", spShort, " on ", island,
                   " are available. Loading results..."))
    terrestrials_Models <- qs::qread(file = terrestrials_models)
    terrestrials_models_time <- qs::qread(file = terrestrials_models_t)
  }
  
  message(paste0("It takes ", round(sum(terrestrials_models_time), 0), 
                 " minutes to run all models."))
  
  modelSelected <- modSel(terrestrials_Models)
  terrestrials_best_model_name <- modelSelected@Full[1,"model"]
  terrestrials_best_model <- terrestrials_Models@fits[[terrestrials_best_model_name]]
  
  ######
  
  ###### MODELS K #####
  
  # We need to check if our K is enough. This is done my running the best model 
  # with various values of K until AIC remains unchanged
  
  terrestrials_best <- paste0("./outputs/",spShort,"_",islandShort,"_bestModel.qs")
  terrestrials_models_K_t <- paste0("./outputs/",spShort,"_",islandShort,"_models_K_time.qs")
  kvals <- c(50, 100, 250, 500, 600, 800) # Min K is obs+1. Crabs have up to 37.
  if (any(reRunK,
          rerunModels, 
          !file.exists(terrestrials_best))){
    if (any(reRunK, rerunModels))
      message(paste0("The argument rerunModels = TRUE or rerunK = TRUE. ", 
                     "Running the models for for assessing K will take some time...")) else   
                       message(paste0("Models for for assessing K for terrestrials on ",
                                      island,
                                      " haven't yet been run yet.",
                                      " Running the models will take some time..."))

    bestLambda <- as.character(mods[models == terrestrials_best_model_name, lambda])
    bestGamma <- as.character(mods[models == terrestrials_best_model_name, gamma])
    bestP <- as.character(mods[models == terrestrials_best_model_name, p])
    bestIota <- if (spShort == "elaenia") as.character(mods[models == terrestrials_best_model_name, iota]) else "iota(.)"
      t0 <- Sys.time()
      en <- environment()
      terrestrials_bestModel <- lapply(1:length(kvals), function(i){
        tic(paste0("Model with K ", kvals[i], " finished: "))
        indMod <- paste0("./outputs/individualModels/", spShort, "_", islandShort, "_K", kvals[i],".qs")
        if (!file.exists(indMod)){
          message(paste0("The model with K ", kvals[i], " doesn't exist. Running."))
        modK <- pcountOpen(lambdaformula = paramTable[variable == bestLambda, value][[1]],  # Initial abundance
                           gammaformula = paramTable[variable == bestGamma, value][[1]],  # Formula for population growth rate
                           omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                           pformula = paramTable[variable == bestP, value][[1]],
                           iotaformula = paramTable[variable == bestIota, value][[1]],
                           data = terrestrialsDF,
                           mixture = modType,
                           dynamics = "trend", # We want the population trend through time
                           immigration = immg,
                           K = kvals[i])
        toc()
        qs::qsave(x = modK, file = indMod)
        } else {
          message(paste0("The model with K ", kvals[i], " already exists. Returning."))
          modK <- qs::qread(indMod)
        }
        assign(x = paste0("t", i), value = Sys.time(), envir = en)
        return(modK)
      })
    names(terrestrials_bestModel) <- paste0("K_", kvals)
    qs::qsave(x = terrestrials_bestModel, file = terrestrials_best)
    terrestrials_models_K_time <- length(terrestrials_bestModel)
    for (m in 0:(terrestrials_models_K_time-1)){
      otm <- get(paste0("t", m+1))-get(paste0("t", m))
      tm <- as.numeric(otm)
      terrestrials_models_K_time[m+1] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    qs::qsave(x = terrestrials_models_K_time, file = terrestrials_models_K_t)
  } else {
    message(paste0("Models for assessing K for terrestrials on ", island,
                   " are available. Loading results..."))
    terrestrials_bestModel <- qs::qread(file = terrestrials_best)
    terrestrials_models_K_time <- qs::qread(file = terrestrials_models_K_t)
  }
  terrestrials_AICs <- sapply(terrestrials_bestModel, function(m) m@AIC)
  
  # Plot AIC to see if the K is ok
  p <- plot(x = kvals, y = terrestrials_AICs, xlab = "K values", ylab = "AIC")

  # We then override the best model with the model with highest K --> most stable parameters
  terrestrials_best_model <- terrestrials_bestModel[[length(terrestrials_bestModel)]]
  ######
  
  ###### PARAMETER ESTIMATION ######
  # Back-transformation of parameters using predict()
  initialPopulation <- predict(obj = terrestrials_best_model, type = "lambda")
  populationGrowthRate <- predict(terrestrials_best_model, type = "gamma")
  detection <- predict(terrestrials_best_model, type = "det")
  immigration <- if (spShort == "elaenia") predict(terrestrials_best_model, type = "iota") else NA
  
  pop <- unique(data.table(initialPop = predict(terrestrials_best_model, type = "lambda")[,1], 
                           cover = terrestrials_best_model@data@siteCovs[, "cover"]))
  
  # The observation radius for terrestrials is 3m, so the area is 28.27m2 or 0,002827ha
  pop[, areaM := (1*pi*area^2)] # height*pi*observation radius^2 = observed area --> 
  # Observations for mabuia were done in 3D, trees, rocks, bushes, etc.
  pop[, densityM := initialPop/areaM]
  ######
  
  return(list(IndividualModelRunTime = terrestrials_models_time,
              totalModelRunTime = sum(terrestrials_models_time),
              nullModels = nullModels,
              allModels = terrestrials_Models,
              bestModelName = terrestrials_best_model_name,
              bestModel = terrestrials_best_model,
              allKmodels = terrestrials_bestModel,
              Kplot = p,
              initialPopulation = initialPopulation,
              detectionProbability = detection,
              populationGrowthRate = populationGrowthRate,
              immigration = immigration,
              initialPopulationDensity = pop))
}
#################################
