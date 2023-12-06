plotPopulationGrowth <- function(wildlifeDataset, ratsDataset, 
                                 whatToPlot = "populationGrowthRate"){
  colsToSelect <- names(ratsDataset)
  DT <- rbindlist(list(unique(wildlifeDataset[variable == whatToPlot, ..colsToSelect]),
                       ratsDataset), use.names = TRUE)
  if (length(whatToPlot) > 2) stop("Only one variable allowed")
  plotDT <- rbindlist(lapply(unique(wildlifeDataset$species), function(sp){
    spDT <- DT[species == sp,]
    availableCycles <- unique(round(spDT$cycle, 0))
    cyclesDT <- rbindlist(lapply(availableCycles, function(cy){
      cyDT <- spDT[round(cycle, 0) == cy,]
      # Estimate
      estimateMeio <- unique(mean(cyDT[island == "Ilha do Meio", estimate]))
      estimateRata <- unique(mean(cyDT[island == "Ilha Rata", estimate]))
      # LCI
      lciRata <- mean(cyDT[island == "Ilha Rata", LCI])-estimateRata
      uciMeio <- mean(cyDT[island == "Ilha do Meio", UCI])-estimateRata
      # UCI
      lciMeio <- mean(cyDT[island == "Ilha do Meio", LCI])-estimateRata
      uciRata <- mean(cyDT[island == "Ilha Rata", UCI])-estimateRata
      
      # Difference islands
      estimateDiff <- mean(estimateMeio-estimateRata)
      estR <- if (all(whatToPlot %in% c("populationGrowthRate",
                                        "immigration"), cy == 1)) NA else 0
      finalDT <- unique(data.table(species = rep(sp, times = 2),
                            cycle = rep(cy, times = 2),
                            island = c("Ilha do Meio", "Ilha Rata"),
                            variable = rep(whatToPlot, times = 2),
                            estimate = c(estimateDiff, estR),
                            UCI = c(uciMeio, uciRata),
                            LCI = c(lciMeio, lciRata)))
      # for (j in seq_len(ncol(finalDT)))
      #   set(finalDT,which(is.na(finalDT[[j]])),j,0) # Not sure I want to replace NA's...
    }))
  }))
  
  # Plotting
  spWd <- sort(unique(wildlifeDataset$species))
  spR <- sort(unique(ratsDataset$species))
  colsToUse <- names(plotDT)
  finalDT <- rbindlist(list(plotDT,
                            ratsDataset[, ..colsToUse]), use.names = TRUE)
  finalDT[, species := factor(species, levels = c(spWd, spR))]
  browser()
  # TRY OTHER MEASURES BECAUSE OF THE FIRST CYCLE OR INSTEAD OF DOING IT 'NA', MAKE THE LAST ONE 'NA'
  # Fix Elaenia. Models are garbage!! Maybe without immigration?
  plotCols <- c("darkgreen", "lightblue")
  p <- ggplot(finalDT, aes(x = cycle, y = estimate)) +
      geom_line(aes(color = island), lwd = 1.5) + 
      scale_color_manual(values = plotCols) +
      scale_fill_manual(values = plotCols) +
      # geom_smooth(aes(color = island), lwd = 1.5, method = "loess") +
      geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = island), alpha = 0.3) +
      geom_point() +
      facet_grid(rows = vars(species), scales = "free_y") +
    theme(legend.position="bottom", legend.title = element_blank())
  p
}
