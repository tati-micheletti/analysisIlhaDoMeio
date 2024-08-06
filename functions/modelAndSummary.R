modelAndSummary <- function(ROW,
                            fullComb,
                            Data,
                            weatherData,
                            reRunK = FALSE,
                            rerunModels = FALSE,
                            rewriteTable = FALSE
                            ){
  gc()
  sp <- as.character(fullComb[ROW, Species])
  isl <- as.character(fullComb[ROW, Island])
  shrt <- as.character(fullComb[ROW, shortSpIsland])
  useRbd <- fullComb[ROW, rbstDsgn]
  useImm <- fullComb[ROW, immgrt]
  Rbd <- if (useRbd) "RbstDsgn" else ""
  Imm <- if (useImm) "Imm" else ""
  assign(shrt, value = invisible(modelWildlife(DT = Data,
                                               island = isl, 
                                               species = sp,
                                               weatherData = weatherData,
                                               useImmigration = useImm,
                                               useRobustDesign = useRbd,
                                               rerunModels = rerunModels, 
                                               reRunK = reRunK)))
  print(shrt)
  changeInPopulation <- calculatePopulationChange(species = sp, 
                                                  island = isl,
                                                  shortName = paste(shrt, Rbd, Imm, sep = "_"),
                                                  envir = environment(),
                                                  RobustDesign = useRbd,
                                                  rewriteTable = rewriteTable)
  rm(shrt)
  gc()
  return(changeInPopulation)
}
