createFullModelTable <- function(overwriteTable = FALSE,
                                 finalTable){
  
  # finalTable is used to filter the correct model (best model with stable parameter values)
  
  fullTableName <- file.path("./outputs", "fullModelsTable.csv")
if (any(!file.exists(fullTableName),
        overwriteTable)){
  fullComb <- data.table(expand.grid(Species = c("Trachylepis atlantica",
                                                 "Elaenia ridleyana",
                                                 "Johngarthia lagostoma",
                                                 "Sula dactylatra"),
                                     Island = c("Ilha do Meio", 
                                                "Ilha Rata")))
  fullComb[, shortSpIsland := c("mabuiaMeio", 
                                "elaeniaMeio",
                                "crabMeio", 
                                "maskedBoobyMeio",
                                "mabuiaRata", 
                                "elaeniaRata",
                                "crabRata", 
                                "maskedBoobyRata")]
  fullComb[, rbstDsgn := c(rep(FALSE, times = 1), 
                           rep(TRUE, times = 1),
                           rep(FALSE, times = 3),
                           rep(TRUE, times = 1),
                           rep(FALSE, times = 2)
  )]
  # fullComb[, immgrt := rep(FALSE, times = 8)]
  
  paramTable <- loadParametersTable() 
  
  allSpMods <- invisible(rbindlist(lapply(1:NROW(fullComb), function(ROW){
    
    Species <- as.character(fullComb[ROW, Species])
    Island <- as.character(fullComb[ROW, Island])
    useRbd <- fullComb[ROW, rbstDsgn]
    # useImm <- fullComb[ROW, immgrt]
    Rbd <- if (useRbd) "RbstDsgn" else ""
    # Imm <- if (useImm) "Immig" else ""
    
    islandShort <-  if (Island == "Ilha Rata") "Rata" else "Meio"
    spShort <-  if (Species == "Elaenia ridleyana") "elaenia" else 
      if (Species == "Trachylepis atlantica") "mabuia" else 
        if (Species == "Johngarthia lagostoma") "crab" else
          if (Species == "Sula dactylatra") "maskedBooby" else
            stop(paste0("Species ", Species, " is not in the dataset"))
    
    # imm <- if (useImm) spShort else NULL
    
    spParameters <- spParams(islandShort = islandShort, 
                             spShort = spShort 
                             # addImmigrationFor = imm
                             )
    
    mods <- data.table(expand.grid(lambda = spParameters$lambda, 
                                   gamma = spParameters$gamma, 
                                   p = spParameters$p, 
                                   iota = spParameters$iota))
    mods[, models := paste0(lambda, gamma, p, iota)]
    mods[, c("lambda", "gamma", "p", "iota") := NULL]
    mods[, Species := spShort]
    mods[, Island := islandShort]
    
    addOns <- if (useRbd) "RbstDsgn" else NULL

    individualModels <- paste("individualModels", addOns, sep = "_")
    
    mods[, modelFileName := file.path("./outputs", individualModels,
                                      paste0(spShort, "_", islandShort, "_", models,".qs"))]
    
    TB <- rbindlist(lapply(1:NROW(mods), function(ROW2){
      
      # Check if this model is in the list of best models. If so, get the correct one! 
      modToLoad <- mods[ROW2, modelFileName]
      modelFormulation <- unique(finalTable[species == Species & island == Island, modelFormulation])
      isBestModel <- if (basename(modToLoad) == paste0(spShort, "_", islandShort, "_", modelFormulation,".qs")) TRUE else FALSE
      if (isBestModel){
        # Update modToLoad to the best models with stable parameter values
        modToLoad <- file.path("./outputs", individualModels,
                               paste0(spShort, "_", islandShort, "_bestModel.qs"))
      }
      modLoaded <- qs::qread(modToLoad)
      if (isBestModel){
        # Need to get only the last one
        if (is(modLoaded, "list"))
          modLoaded <- modLoaded[[length(modLoaded)]]
      }
      mdSumm <- suppressMessages(invisible(summary(modLoaded)))
      allPars <- names(mdSumm)
      
      pars <- rbindlist(lapply(allPars, function(Ps){
        paramName <- Ps
        covariate <- rownames(mdSumm[[Ps]])
        estimate <- mdSumm[[Ps]][["Estimate"]]
        SE <- mdSumm[[Ps]][["SE"]]
        z <- mdSumm[[Ps]][["z"]]
        p <- mdSumm[[Ps]][["P(>|z|)"]]
        return(data.table(parameter = paramName,
                          covariate = covariate,
                          estimate = estimate,
                          SE = SE,
                          z = z,
                          p = p))
      }))
      
      tbDT <- data.table(species = Species,
                         island = Island,
                         model = mods[ROW2, models],
                         parameter = pars[["parameter"]],
                         covariate = pars[["covariate"]],
                         estimate = pars[["estimate"]],
                         standardError = pars[["SE"]],
                         z = pars[["z"]],
                         p = pars[["p"]],
                         AIC = modLoaded@AIC,
                         modelMixture = modLoaded@mixture,
                         K = modLoaded@K
                         # immigration = modLoaded@immigration
      )
      return(tbDT)
    }))
    
    return(TB)
    
  })))
  
  write.csv(x = allSpMods, file = fullTableName)
  
} else {
  
  allSpMods <- data.table::fread(fullTableName)

  }

return(allSpMods)
}

