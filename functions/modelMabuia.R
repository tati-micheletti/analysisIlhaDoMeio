################################# MABUIA #################################

modelMabuia <- function(DT, island, rerunModels = FALSE){

  # First we subset the data to the interested place and the species
Mabuia <- DT[Island == island & Species == "Trachylepis atlantica",]
area <- unique(Mabuia$Radius)
# Drop what we don't need, or don't want: 
# Island, Year, Month, Day, Species, Radius, Original_Landscape, originDate
toRemove <- c("Island", "Year", "Month", "Day", "Species", 
              "Radius", "Original_Landscape", "originDate",
              "pointID", "transectID") # For mabuia in Rata, we need to use ponto_numero as site
                                       # because there were mistakes in the table that were corrected
                                       # by PM later.
if (island == "Ilha Rata"){
  toRemove <- c(toRemove, "site")
  Mabuia <- Mabuia[, (toRemove) := NULL]
  names(Mabuia)[names(Mabuia)=="ponto_numero"] <- "site"
  Mabuia[, site := as.character(site)]
} else {
  toRemove <- c(toRemove, "ponto_numero")
  Mabuia <- Mabuia[, (toRemove) := NULL]
}

# Now we make the needed objects:
# y --> Matrix with observations where rows represent each site (unmarked calls this 'M'), and columns the 
# (increasing) sampling occasion (i.e., the time series of counts, representing 
# each visit), which unmarked calls 'T'. The value is the Counts.
finalCounts <- dcast(Mabuia, site ~ Date, value.var = "Counts")
y <- as.matrix(finalCounts[,-1]) # Convert to matrix

# obsCov --> Observation covariates are the ones that vary per sampling occasion.
# This object is a named list of each covariate, which are data frames of rows 
# representing each site and the columns each sampling occasion. 
# Examples are: observerID, JulianDate, timeSinceStartEradication
observerID <- as.matrix(dcast(Mabuia, site ~ Date, value.var = "Observer")[,-1])
JulianDate <- as.matrix(dcast(Mabuia, site ~ Date, value.var = "JulianDate")[,-1])
timeSinceStartEradication <- as.matrix(dcast(Mabuia, site ~ Date, value.var = "timeSinceStartEradication")[,-1])

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

obsCovs <- list(observerID = observerID,
                JulianDate = JulianDate,
                timeSinceStartEradication = timeSinceStartEradication)

# siteCovs --> Site covariates are fixed through time, varying only among each 
# other. This object is matrix where rows representing each site and the columns 
# each covariate. 
# Examples are: Cover (landscape type)
cover <- data.frame(dcast(Mabuia, site ~ Date, value.var = "Landscape")[,4]) 
#TODO Here I need to make sure all rows have information and then I can get whatever row
names(cover) <- "cover"

# numPrimary --> Number of primary periods
numPrimary <- NCOL(y)

yearlySiteCovs <- list(JulianDate = JulianDate,
                       timeSinceStartEradication = timeSinceStartEradication) 

# Create the unmarkedFrame
mabuiaDF <- unmarkedFramePCO(y = y,
                             siteCovs = cover,
                             obsCovs = obsCovs,
                             numPrimary = numPrimary,
                             yearlySiteCovs = yearlySiteCovs
)
# Take a look
mabuiaDF
summary(mabuiaDF)

# MODELS TESTED
# Null Model
# Here we test for Zero Inflation on the data
test_ZI(dts = na.omit(as.numeric(y)))

# We run then a negative binomial (NB) and a zero inflated poisson (ZIP)
# We don't run Poisson because the data us zero inflated
if (island == "Ilha do Meio"){
  mabuia_nullModels_file <- "./outputs/mabuia_Meio_nullModels.qs"  
} else {
  mabuia_nullModels_file <- "./outputs/mabuia_Rata_nullModels.qs"  
}

if (any(rerunModels, !file.exists(mabuia_nullModels_file))){
  if (rerunModels)
    message(paste0("The argument rerunModels = TRUE.", 
                   "Running the null models should not take too much time...")) else   
                     message(paste0("Null models for Mabuia on ", island," haven't yet been run yet.",
                                    " Running these models should not take too much time..."))
t0 <- Sys.time()
nullNB <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                     gammaformula = ~1,  # Formula for population growth rate
                     omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                     pformula = ~1,
                     data = mabuiaDF,
                     mixture = "NB", # Negative binomial
                     dynamics = "trend", # We want the population trend through time
                     immigration = FALSE,
                     K = 100) # AIC: 648.5343
# The dispersion of data is not signifficant

nullZIP <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                      gammaformula = ~1,  # Formula for population growth rate
                      omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                      pformula = ~1,
                      data = mabuiaDF,
                      mixture = "ZIP", # Negative binomial
                      dynamics = "trend", # We want the population trend through time
                      immigration = FALSE,
                      K = 100) # AIC: 707.3251
# We do have Zero Inflated data (significant)
nullModels <- list(nullZIP = nullZIP,
                   nullNB = nullNB)
qs::qsave(x = nullModels, file = mabuia_nullModels_file)
} else {
  message(paste0("Null models for Mabuia on ", island," are available. Loading results..."))
  nullModels <- qs::qread(file = mabuia_nullModels_file)
}
# This points us to use ZIP models instead of a NB 
# Supported by Stoklosa et al., (Diversity 2022, 14(5), 320; https://doi.org/10.3390/d14050320) 
# and by Kery 2018 (Ecology, 99(2), 2018, pp. 281–288), even though it has smaller AIC for NB
if (island == "Ilha do Meio"){
  mabuia_models <- "./outputs/mabuia_Meio_models.qs"
  mabuia_models_t <- "./outputs/mabuia_Meio_models_time.qs"
} else {
  mabuia_models <- "./outputs/mabuia_Rata_models.qs"
  mabuia_models_t <- "./outputs/mabuia_Rata_models_time.qs"
}
if (any(rerunModels, !file.exists(mabuia_models))){
  if (rerunModels)
    message(paste0("The argument rerunModels = TRUE.", 
                   "Running the models will take some time...")) else   
      message(paste0("Models for Mabuia on ", island," haven't yet been run yet.",
                     " Running the models will take some time..."))
t1 <- Sys.time()
mod1 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~1,  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~1,
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100)
t2 <- Sys.time()
mod2 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                   gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~1,
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100)
t3 <- Sys.time()
mod3 <- pcountOpen(lambdaformula = ~1,  # Initial abundance
                   gammaformula = ~1,  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~scale(timeSinceStartEradication),
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t4 <- Sys.time()
mod4 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~1,
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t5 <- Sys.time()
mod5 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~1,  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~scale(timeSinceStartEradication),
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t6 <- Sys.time()
mod6 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~scale(timeSinceStartEradication),
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t7 <- Sys.time()
mod7 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~scale(JulianDate),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~1,
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t8 <- Sys.time()
mod8 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~scale(timeSinceStartEradication)+scale(JulianDate),
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t9 <- Sys.time()
mod9 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                   gammaformula = ~scale(timeSinceStartEradication)+scale(JulianDate),  # Formula for population growth rate
                   omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                   pformula = ~scale(timeSinceStartEradication)+scale(JulianDate)+observerID,
                   data = mabuiaDF,
                   mixture = "ZIP", # Negative binomial
                   dynamics = "trend", # We want the population trend through time
                   immigration = FALSE,
                   K = 100) 
t10 <- Sys.time()
mod10 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(timeSinceStartEradication)+scale(JulianDate),
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t11 <- Sys.time()
mod11 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(timeSinceStartEradication)+scale(JulianDate)+observerID,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t12 <- Sys.time()
mod12 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(JulianDate)+observerID,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 

t13 <-Sys.time()
mod13 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(JulianDate),
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t14 <- Sys.time()
mod14 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(JulianDate)+cover,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t15 <-Sys.time()
mod15 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~cover,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t16 <- Sys.time()
mod16 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~cover+scale(timeSinceStartEradication),  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(JulianDate)+observerID,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t17 <- Sys.time()
mod17 <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                    gammaformula = ~1,  # Formula for population growth rate
                    omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                    pformula = ~scale(JulianDate)+observerID,
                    data = mabuiaDF,
                    mixture = "ZIP", # Negative binomial
                    dynamics = "trend", # We want the population trend through time
                    immigration = FALSE,
                    K = 100) 
t18 <-Sys.time()

mabuia_Models <- fitList(
  "lam(.)gamma(.)phi(.)p(.)" = nullZIP,
  "lam(cover)gamma(.)phi(.)p(.)" = mod1,
  "lam(.)gamma(timeSinceStartEradication)phi(.)p(.)" = mod2,  
  "lam(.)gamma(.)phi(.)p(timeSinceStartEradication)" = mod3, 
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(.)" = mod4, 
  "lam(cover)gamma(.)phi(.)p(timeSinceStartEradication)" = mod5,  
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(timeSinceStartEradication)" = mod6,  
  "lam(cover)gamma(JulianDate)phi(.)p(.)" = mod7, 
  "lam(cover)gamma(timeSinceStartEradication+JulianDate)phi(.)p(timeSinceStartEradication+JulianDate)" = mod8, 
  "lam(cover)gamma(timeSinceStartEradication+JulianDate)phi(.)p(timeSinceStartEradication+JulianDate+observer)" = mod9,
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(timeSinceStartEradication+JulianDate)" = mod10,
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(timeSinceStartEradication+JulianDate+observer)" = mod11,
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(JulianDate+observer)" = mod12, # BEST MODEL
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(JulianDate)" = mod13,
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(JulianDate+cover)" = mod14,
  "lam(cover)gamma(timeSinceStartEradication)phi(.)p(cover)" = mod15,
  "lam(cover)gamma(cover+timeSinceStartEradication)phi(.)p(JulianDate+observer)" = mod16,
  "lam(cover)gamma(.)phi(.)p(JulianDate+observer)" = mod17
)

mabuia_models_time <- numeric(length(mabuia_Models@fits))
for (m in 0:(length(mabuia_models_time)-1)){
  otm <- get(paste0("t", m+1))-get(paste0("t", m))
  tm <- as.numeric(otm)
  mabuia_models_time[m+1] <- if (units(otm)=="secs")
    tm/60 else # Seconds
      if (units(otm)=="mins")
        tm else # Minutes
          tm*60 # Hours
}

qs::qsave(x = mabuia_Models, file = mabuia_models)
qs::qsave(x = mabuia_models_time, file = mabuia_models_t)

} else {
  message(paste0("Models for Mabuia on ", island," are available. Loading results..."))
  mabuia_Models <- qs::qread(file = mabuia_models)
  mabuia_models_time <- qs::qread(file = mabuia_models_t)
}

message(paste0("It takes ", round(sum(mabuia_models_time), 0), 
               " minutes to run all models."))

modelSelected <- modSel(mabuia_Models)
mabuia_best_model_name <- modelSelected@Full[1,"model"]
mabuia_best_model <- mabuia_Models@fits[[mabuia_best_model_name]]

# I won't deal with calculating the total abundance, not the proposal here! We want to see the change in 
# lambda / population growth rate

# We need to check if our K is enough. This is done my running the best model with various values of 
# K until AIC remains unchanged
if (island == "Ilha do Meio"){
  mabuia_best <- "./outputs/mabuia_Meio_bestModel.qs"
} else {
  mabuia_best <- "./outputs/mabuia_Rata_bestModel.qs"  
}
kvals <- c(20, 50, 100, 250, 500, 600, 800) # Min K is 13+1 (max observed number of indiv + 1) so we go between 20 and 500
if (any(rerunModels, 
        !file.exists(mabuia_best))){
  if (rerunModels)
    message(paste0("The argument rerunModels = TRUE.", 
                   "Running the models for for assessing K will take some time...")) else   
                     message(paste0("Models for for assessing K for Mabuia on ", island,
                                    " haven't yet been run yet. ",
                                    "Running the models will take some time..."))
  mabuia_bestModel <- list()
  for (i in 1:length(kvals)){
    tic(paste0("Model with K ", kvals[i], " finished: "))
    if (island == "Ilha do Meio"){
    mabuia_bestModel[[i]] <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                                             gammaformula = ~scale(timeSinceStartEradication),  # Formula for population growth rate
                                             omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                                             pformula = ~scale(JulianDate)+observerID,
                                             data = mabuiaDF2,
                                             mixture = "ZIP", # Negative binomial
                                             dynamics = "trend", # We want the population trend through time
                                             immigration = FALSE,
                                             K = kvals[i])
    } else {
      mabuia_bestModel[[i]] <- pcountOpen(lambdaformula = ~cover,  # Initial abundance
                                          gammaformula = ~1,  # Formula for population growth rate
                                          omegaformula = ~1, # Formula for apparent survival probability: can't use covs here
                                          pformula = ~scale(JulianDate)+observerID,
                                          data = mabuiaDF2,
                                          mixture = "ZIP", # Negative binomial
                                          dynamics = "trend", # We want the population trend through time
                                          immigration = FALSE,
                                          K = kvals[i])
    }
    toc()
  }
  names(mabuia_bestModel) <- paste0("K_", kvals)
  qs::qsave(x = mabuia_bestModel, file = mabuia_best)
} else {
  message(paste0("Models for assessing K for Mabuia on ", island,
                 " are available. Loading results..."))
  mabuia_bestModel <- qs::qread(file = mabuia_best)
}
mabuia_AICs <- sapply(mabuia_bestModel, function(m) m@AIC)

# Plot AIC to see if the K is ok
p <- plot(x = kvals, y = mabuia_AICs, xlab = "K values", ylab = "AIC")
# axis(1, at=mabuia_AICs, labels=names(mabuia_AICs))

# We then override the best model with the model with highest K --> most stable parameters
# mabuia_best_model <- mabuia_bestModel[[length(mabuia_bestModel)]]
mabuia_best_model <- mabuia_bestModel[[length(mabuia_bestModel)]]
# Back-transformation of parameters using predict()
populationGrowthRate <- predict(mabuia_best_model, type = "gamma")
initialPopulation <- predict(obj = mabuia_best_model, type = "lambda")
detection <- predict(mabuia_best_model, type = "det")

# Gamma (lambda --> population growth rate is increasing)
pop <- unique(data.table(initialPop = predict(mabuia_best_model, type = "lambda")[,1], 
                         cover = mabuia_best_model@data@siteCovs[, "cover"]))

# The observation radius for Mabuia is 3m, so the area is 28.27m2 or 0,002827ha
pop[, areaM3 := (1*pi*area^2)] # height*pi*observation radius^2 = observed area --> Observations were done in 3D, trees, rocks, bushes, etc.
pop[, densityM3 := initialPop/areaM3]
# 0.357 ± 0.170 individuals/m² (Vini's work) on secondary islands (0.187 - 0.527)

return(list(IndividualModelRunTime = mabuia_models_time,
            totalModelRunTime = sum(mabuia_models_time),
            nullModels = nullModels,
            allModels = mabuia_Models,
            bestModelName = mabuia_best_model_name,
            bestModel = mabuia_best_model,
            allKmodels = mabuia_bestModel,
            Kplot = p,
            initialPopulation = initialPopulation,
            detectionProbability = detection,
            populationGrowthRate = populationGrowthRate,
            initialPopulationDensity = pop))
}
#################################
