################################# ELAENIA ###########################

modelElaenia <- function(DT, island, rerunModels = FALSE, reRunK = FALSE){
  # First we subset the data to the interested place and the species
  elenia <- DT[Island == island & Species == "Elaenia ridleyana",]
  area <- unique(elenia$Radius)
  # Drop what we don't need, or don't want: 
  # Island, Year, Month, Day, Species, Radius, Original_Landscape, originDate
  toRemove <- c("Island", "Year", "Month", "Day", "Species", 
                "Radius", "Original_Landscape", "originDate",
                "pointID", "transectID")
  elenia <- elenia[, (toRemove) := NULL]
  
  # Now we make the needed objects:
  # y --> Matrix with observations where rows represent each site (unmarked calls this 'M'), and columns the 
  # (increasing) sampling occasion (i.e., the time series of counts, representing 
  # each visit), which unmarked calls 'T'. The value is the Counts.
  finalCounts <- dcast(elenia, site ~ Date, value.var = "Counts")
  y <- as.matrix(finalCounts[,-1]) # Convert to matrix
  
  # obsCov --> Observation covariates are the ones that vary per sampling occasion.
  # This object is a named list of each covariate, which are data frames of rows 
  # representing each site and the columns each sampling occasion. 
  # Examples are: JulianDate, timeSinceStartEradication. # Observer is a moot point for Elenia, 
  # as only one person did it
  JulianDate <- as.matrix(dcast(elenia, site ~ Date, value.var = "JulianDate")[,-1])
  timeSinceStartEradication <- as.matrix(dcast(elenia, site ~ Date, value.var = "timeSinceStartEradication")[,-1])
  
  # Fill in the matrix of JulianDate and timeSinceStartEradication: this affects the models later on
  replaceNAwithColMean  <- function(mx){
    DT <- as.data.table(mx)
    nm <- names(DT)[colSums(is.na(DT)) != 0]
    for(clnm in nm){
      set(DT, i = which(is.na(DT[[clnm]])), j = clnm, 
          value = mean(DT[, get(clnm)], na.rm = TRUE))
    }
    return(as.matrix(DT))
  }
  JulianDate <- replaceNAwithColMean(JulianDate)
  timeSinceStartEradication <- replaceNAwithColMean(timeSinceStartEradication)
  
  obsCovs <- list(JulianDate = JulianDate,
                  timeSinceStartEradication = timeSinceStartEradication)
  
  # siteCovs --> Site covariates are fixed through time, varying only among each 
  # other. This object is matrix where rows representing each site and the columns 
  # each covariate. 
  # Examples are: Cover (landscape type)
  cover <- data.frame(dcast(elenia, site ~ Date, value.var = "Landscape")[,2])
  names(cover) <- "cover"
  
  numPrimary <- NCOL(y)/2 # In the case of Elenia, we have 2 counts (i.e., secondary 
  # periods) per primary occasion
  
  # Since the closed periods are really close, we can use one variable for both days 
  # (yearlySiteCovs). To do this, we need to get the mean of every n columns (n = 
  # secondary periods, or 2 in our case). I will do this using a loop
  interval <- NCOL(y)/numPrimary # How many secondary periods we have
  firstColsToJoin <- seq(from = 1, to = NCOL(y), by = interval)
  # For timeSinceStartEradication
  newTSE <- lapply(firstColsToJoin, function(i){
    colsToJoin <- colnames(timeSinceStartEradication)[i:(i+1)]
    TSE <- data.table(timeSinceStartEradication[, colsToJoin])
    TSE[, colnames(timeSinceStartEradication)[i] := rowMeans(.SD), .SDcols = names(TSE)]
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
                         timeSinceStartEradication = as.matrix(newTSE)) 
  
  # Create the unmarkedFrame
  eleniaDF <- unmarkedFramePCO(y = y,
                               siteCovs = cover,
                               obsCovs = obsCovs,
                               numPrimary = numPrimary,
                               yearlySiteCovs = yearlySiteCovs
  )
  # Take a look
  eleniaDF
  summary(eleniaDF)
  
  # MODELS TESTED
  # Null Model
  # Here we test for Zero Inflation on the data
  test_ZI(dts = na.omit(as.numeric(y)))
  
  # We run then a negative binomial (NB) and a Poisson (ZIP)
  # We don't run ZIP because the data is NOT zero inflated
  if (island == "Ilha do Meio"){
    elenia_nullModels_file <- "./outputs/elenia_Meio_nullModels.qs"  
  } else {
    elenia_nullModels_file <- "./outputs/elenia_Rata_nullModels.qs"  
  }
  
  if (any(rerunModels, !file.exists(elenia_nullModels_file))){
    if (rerunModels)
      message(paste0("The argument rerunModels = TRUE.", 
                     "Running the null models should not take too much time...")) else   
                       message(paste0("Null models for Elaenia on ", island, " haven't yet been run yet.",
                                      " Running these models should not take too much time..."))
    t0 <- Sys.time()
    nullP <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                        gammaformula = ~1,  # Formula for population growth rate
                        omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                        pformula = ~1,
                        data = eleniaDF,
                        mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                        dynamics = "trend", # We want the population trend through time
                        iotaformula = ~1,
                        immigration = TRUE,
                        K = 100) # AIC: 227.2522
    
    nullNB <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                         gammaformula = ~1,  # Formula for population growth rate
                         omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                         pformula = ~1,
                         data = eleniaDF,
                         iotaformula = ~1,
                         mixture = "NB", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                         dynamics = "trend", # We want the population trend through time
                         immigration = TRUE,
                         K = 100) # AIC: 229.2535 
    
    # We do have Zero Inflated data (significant)
    nullModels <- list(nullP = nullP,
                       nullNB = nullNB)
    qs::qsave(x = nullModels, file = elenia_nullModels_file)
  } else {
    message(paste0("Null models for Elaenia on ", island, 
                   " are available. Loading results..."))
    nullModels <- qs::qread(file = elenia_nullModels_file)
  }
  # This points us to use Poisson models
  if (island == "Ilha do Meio"){
    elenia_models <- "./outputs/elenia_Meio_models.qs"
    elenia_models_t <- "./outputs/elenia_Meio_models_time.qs"
  } else {
    elenia_models <- "./outputs/elenia_Rata_models.qs"
    elenia_models_t <- "./outputs/elenia_Rata_models_time.qs"
  }
  if (any(rerunModels, !file.exists(elenia_models))){
    if (rerunModels)
      message(paste0("The argument rerunModels = TRUE.", 
                     "Running the models will take some time...")) else   
                       message(paste0("Models for Elaenia ", island, " haven't yet been run yet.",
                                      " Running the models will take some time..."))
    t1 <- Sys.time()
    mod1 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       data = eleniaDF,
                       iotaformula = ~1,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t2 <- Sys.time()
    mod2 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~1,  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~scale(timeSinceStartEradication),
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t3 <- Sys.time()
    mod3 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~scale(timeSinceStartEradication),
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t4 <- Sys.time()
    mod4 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(JulianDate),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~1,
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t5 <- Sys.time()
    mod5 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~1,  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~scale(JulianDate),
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t6 <- Sys.time()
    mod6 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(JulianDate),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~scale(JulianDate),
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t7 <- Sys.time()
    mod7 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~1,
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t8 <- Sys.time()
    mod8 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~1,  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       data = eleniaDF,
                       iotaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:
    t9 <- Sys.time()
    mod9 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                       gammaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),  # Formula for population growth rate
                       omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                       pformula = ~1,
                       iotaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),
                       data = eleniaDF,
                       mixture = "P", # Not overdispersed and not zero inflated. Here Poisson is the best model.
                       dynamics = "trend", # We want the population trend through time
                       immigration = TRUE,
                       K = 100) # AIC:

    t10 <- Sys.time()
    elenia_Models <- fitList(
      "lam(.)gamma(.)phi(.)p(.)iota(.)" = nullP,
      "lam(.)gamma(timeSinceStartEradication)phi(.)p(.)iota(.)" = mod1,
      "lam(.)gamma(.)phi(.)p(.)iota(timeSinceStartEradication)" = mod2,  # Best model
      "lam(.)gamma(timeSinceStartEradication)phi(.)p(.)iota(timeSinceStartEradication)" = mod3,  # 2nd Best model
      "lam(.)gamma(JulianDate)phi(.)p(.)iota(.)" = mod4,
      "lam(.)gamma()phi(.)p(.)iota(JulianDate)" = mod5,
      "lam(.)gamma(JulianDate)phi(.)p(.)iota(JulianDate)" = mod6,
      "lam(.)gamma(timeSinceStartEradication+JulianDate)phi(.)p(.)iota(.)" = mod7,
      "lam(.)gamma(.)phi(.)p(.)iota(timeSinceStartEradication+JulianDate)" = mod8,
      "lam(.)gamma(timeSinceStartEradication+JulianDate)phi(.)p(.)iota(timeSinceStartEradication+JulianDate)" = mod9
    )

    elenia_models_time <- numeric(length(elenia_Models@fits))
    for (m in 0:(length(elenia_models_time)-1)){
      otm <- get(paste0("t", m+1))-get(paste0("t", m))
      tm <- as.numeric(otm)
      elenia_models_time[m+1] <- if (units(otm)=="secs")
        tm/60 else # Seconds
          if (units(otm)=="mins")
            tm else # Minutes
              tm*60 # Hours
    }
    qs::qsave(x = elenia_Models, file = elenia_models)
    qs::qsave(x = elenia_models_time, file = elenia_models_t)
  } else {
    message(paste0("Models for Elaenia on ", island, 
                   " are available. Loading results..."))
    elenia_Models <- qs::qread(file = elenia_models)
    elenia_models_time <- qs::qread(file = elenia_models_t)
  }
  
  message("It takes ", round(sum(elenia_models_time), 0), " minutes to run all models.")

  modelSelected <- modSel(elenia_Models)
  Elaenia_best_model_name <- modelSelected@Full[1,"model"]
  Elaenia_best_model <- elenia_Models@fits[[Elaenia_best_model_name]]

  # I won't deal with calculating the total abundance, not the proposal here! We want to see the change in 
  # lambda / population growth rate
  
  # We need to check if our K is enough. This is done my running the best model with various values of 
  # K until AIC remains unchanged
  if (island == "Ilha do Meio"){
    elenia_best <- "./outputs/elenia_Meio_bestModel.qs"
  } else {
    elenia_best <- "./outputs/elenia_Rata_bestModel.qs"
  }
  kvals <- c(20, 50, 100, 250, 500, 600, 800) # Min K is 13+1 (max observed number of indiv + 1) so we go between 20 and 500
  if (any(reRunK,
          rerunModels, 
          !file.exists(elenia_best))){
    if (any(reRunK, rerunModels))
      message(paste0("The argument rerunModels = TRUE or rerunK = TRUE. ", 
                     "Running the models for for assessing K will take some time...")) else   
                       message(paste0("Models for for assessing K for Elaenia on Meio Island",
                                      " haven't yet been run yet.",
                                      " Running the models will take some time..."))
    elenia_bestModel <- list()
    for (i in 1:length(kvals)){
      tic(paste0("Model with K ", kvals[i], " for ", island, " finished: "))
      if (island == "Ilha do Meio"){
        elenia_bestModel[[i]] <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                                            gammaformula = ~1,  # Formula for population growth rate
                                            omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                                            pformula = ~1,
                                            iotaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),
                                            data = eleniaDF,
                                            starts = coef(Elaenia_best_model), # To help converge quicker, pass the best model's coefficients
                                            mixture = "P", # Negative binomial
                                            dynamics = "trend", # We want the population trend through time
                                            immigration = TRUE,
                                            K = kvals[i])
      } else {
        elenia_bestModel[[i]] <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                                            gammaformula = ~1,  # Formula for population growth rate
                                            omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                                            pformula = ~1,
                                            iotaformula = ~scale(JulianDate),
                                            starts = coef(Elaenia_best_model),
                                            data = eleniaDF,
                                            mixture = "P", # Negative binomial
                                            dynamics = "trend", # We want the population trend through time
                                            immigration = TRUE,
                                            K = kvals[i])
      }
      toc()
    }
    names(elenia_bestModel) <- paste0("K_", kvals)
    qs::qsave(x = elenia_bestModel, file = elenia_best)
  } else {
    message(paste0("Models for assessing K for Elaenia on ", island,
                   " are available. Loading results..."))
    elenia_bestModel <- qs::qread(file = elenia_best)
  }
  elenia_AICs <- sapply(elenia_bestModel, function(m) m@AIC)
  
  # Plot AIC to see if the K is ok
  p <- plot(x = kvals, y = elenia_AICs, xlab = "K values", ylab = "AIC")
  # axis(1, at=elenia_AICs, labels=names(elenia_AICs))
  
  # We then override the best model with the model with highest K --> most stable parameters
  # Elaenia_best_model <- elenia_bestModel[[length(elenia_bestModel)]]
  Elaenia_best_model <- elenia_bestModel[[length(elenia_bestModel)]]
  # Back-transformation of parameters using predict()
  populationGrowthRate <- predict(Elaenia_best_model, type = "gamma")
  initialPopulation <- predict(obj = Elaenia_best_model, type = "lambda")
  detection <- predict(Elaenia_best_model, type = "det")
  immigration <- predict(Elaenia_best_model, type = "iota")
  
  # Gamma (lambda --> population growth rate is increasing)
  pop <- unique(data.table(initialPop = predict(Elaenia_best_model, type = "lambda")[,1], 
                           cover = Elaenia_best_model@data@siteCovs[, "cover"]))
  
  # The observation radius for Elaenia is 3m, so the area is 28.27m2 or 0,002827ha
  pop[, areaM2 := (pi*area^2)] # height*pi*observation radius^2 = observed area --> Observations were done in 3D, trees, rocks, bushes, etc.
  pop[, densityM2 := initialPop/areaM2]
  # 0.357 ± 0.170 individuals/m² (Vini's work) on secondary islands (0.187 - 0.527)
  
  return(list(IndividualModelRunTime = elenia_models_time,
              totalModelRunTime = sum(elenia_models_time),
              nullModels = nullModels,
              allModels = elenia_Models,
              bestModelName = Elaenia_best_model_name,
              bestModel = Elaenia_best_model,
              allKmodels = elenia_bestModel,
              Kplot = p,
              initialPopulation = initialPopulation,
              detectionProbability = detection,
              populationGrowthRate = populationGrowthRate,
              initialPopulationDensity = pop))
}

#################################
