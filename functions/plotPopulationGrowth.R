plotPopulationGrowth <- function(wildlifeDataset, ratsDataset, 
                                 whatToPlot = "populationGrowthRate"){
  colsToSelect <- names(ratsDataset)
  DT <- rbindlist(list(unique(wildlifeDataset[parameter == whatToPlot, ..colsToSelect]),
                       ratsDataset), use.names = TRUE)
  if (length(whatToPlot) > 2) stop("Only one parameter allowed")
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
                            parameter = rep(whatToPlot, times = 2),
                            estimate = c(estimateDiff, estR),
                            UCI = c(uciMeio, uciRata),
                            LCI = c(lciMeio, lciRata)))
    }))
  }))
  
  # Plotting
  spWd <- sort(unique(wildlifeDataset$species))
  spR <- sort(unique(ratsDataset$species))
  colsToUse <- names(plotDT)
  finalDT <- rbindlist(list(plotDT,
                            ratsDataset[, ..colsToUse]), use.names = TRUE)
  finalDT[, species := factor(species, levels = c(spWd, spR))]
  # Remove Ilha Rata, as this is plot is about the difference
  finalDT <- finalDT[island == "Ilha do Meio",]
  # Remove NA
  finalDT <- finalDT[!is.na(estimate),]
  # Remove unecessary columns: island, parameter
  finalDT[,c("island", "parameter") := NULL]
  finalDT[, speciesType := fifelse(species == "Rattus rattus", "Invasive Species", "Endangered Species")]
  finalDT[species == "Elaenia ridleyana", cycle := scales::rescale(cycle, to = c(1, 4))]
  # Rename for better plotting
  finalDT[species == "Elaenia ridleyana", species := "E. ridleyana"]
  finalDT[species == "Johngarthia lagostoma", species := "J. lagostoma"]
  finalDT[species == "Sula dactylatra", species := "S. dactylatra"]
  finalDT[species == "Trachylepis atlantica", species := "T. atlantica"]
  finalDT[species == "Rattus rattus", species := "R. rattus"]
  finalDT[, species := factor(species, levels = c("S. dactylatra",
                                            "J. lagostoma",
                                            "E. ridleyana",
                                            "T. atlantica",
                                            "R. rattus"))]
  scales_y <- list(
    "E. ridleyana" = scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.4)),
    "J. lagostoma" = scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.4)),
    "S. dactylatra" = scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.4)),
    "T. atlantica" = scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.4)),
    "R. rattus" = scale_y_continuous(limits = c(0, 1073), breaks = seq(0, 1073, 250))
  )
  plotCols <- c("darkgreen", "blue", "red", "orange", "purple")
  plotCols2 <- rep("black", times = 5)
  # browser()
  p <- ggplot(finalDT, aes(x = cycle, y = estimate, group = species)) + #, linetype = species
      geom_point() + #aes(shape = species)
      geom_smooth(method = "lm",  fullrange = TRUE, se = F, aes(color = species)) +
      xlim(1, 4) +
      scale_color_manual(values = plotCols2) +
      geom_hline(yintercept = 0, linetype = "longdash", color = "darkgrey") +
      scale_fill_manual(values = plotCols2) +
      geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = species), alpha = 0.3) +
      geom_point() +
      facet_grid_sc(rows = vars(species), scales = list(y = scales_y)) +
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    xlab("Cycle") + ylab("Population Density                                                                      Differential Population Growth Index                              ")
  p
  nm <- file.path("outputs/Figure4.png")
  ggsave(device = "png", filename = nm, 
         width = 6, height = 10)
}
