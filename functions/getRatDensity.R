getRatDensity <- function(plotSamplingAreas = FALSE){
  #Before we model rat abundance, we prepare a mask of the island so the analysis
  #occurs only on the island area. This is recommended for island environments due
  #to potential edge effects.
  
  PARNA <- terra::vect(file.path(wd, "Geomorfology.shp")) 
  meio <- subset(PARNA, PARNA$AREA_HA == 17.034712, )
  aug17 <- secr::read.capthist(captfile =  file.path(wd, "CAPTURE1.txt"),
                               trapfile = file.path(wd, "TRAPS1.txt"))
  trapPoints <- traps(aug17)
  trapPointsOK <- vect(as.data.frame(trapPoints), geom=c("x", "y"), crs=crs(meio), keepgeom=FALSE)
  if (plotSamplingAreas){
    plot(main = "First Sampling", meio, col = "lightblue")
    plot(trapPointsOK, add = TRUE, col = "red")
  }
  # We then check if the traps are located correctly in the shapefile.
  clippedmask <- make.mask(trapPoints, type = 'trapbuffer', 
                           buffer = 200, poly = meio)
  par(mfrow = c(1,1), mar = c(1,1,1,1))
  plot(clippedmask, border = 100, ppoly = FALSE)
  gm <- as.data.frame(geom(meio))
  polygon(data.frame(x = gm$x, y = gm$y), col = 'lightgreen', border = NA)
  if (plotSamplingAreas){
    plot(clippedmask, dots = FALSE, mesh = grey(0.4), col = NA, polycol = 'blue', add = TRUE)
    plot(trapPoints, detpar = list(pch = 16, cex = 0.8), add = TRUE)
  }
  
  # Everything seems to be fine. Now we can perform our analysis.
  ## Meio, August 2017 - First monitoring campaign
  invisible({
    aug17 <- secr::read.capthist(captfile =  file.path(wd, "CAPTURE1.txt"),
                                 trapfile = file.path(wd, "TRAPS1.txt")) # We have already read the file before so we could extract the trap's information. Here it is again for consistency.
    g0.1 <- secr.fit(aug17, model = list(g0 ~ 1, sigma ~ 1), CL = TRUE,
                     verify = TRUE, mask = clippedmask)
    
    densityAUG17 <- derived(g0.1)
  })
  # We tested different model formulations (i.e., landscape covariates and also 
  # different buffers; see previous commits to see code).  
  # Model diagnostics (i.e., AIC) supported the use of the simplest model (i.e., 
  # detection probability doesn't change), and the island shapefile as a mask. 
  # Applying the mask or buffer of 200 didn't change the results, but other buffers 
  # did.
  
  #Home range results (sigma estimate, lcl = lowe95%CI; ucl = upper95%CI):
  #summary(g0.1)$predicted
  
  #Density results (D estimate, lcl = lowe95%CI; ucl = upper95%CI):
  #densityAUG17
  
  ## Meio, October 2017 - Second monitoring campaign
  invisible({
    oct17 <- secr::read.capthist(captfile =  file.path(wd, "CAPTURE2.txt"),
                                 trapfile = file.path(wd, "TRAPS2.txt"))
    g0.2 <- secr.fit(oct17, model = list(g0 ~ 1, sigma ~ 1), CL = TRUE,
                     verify = TRUE, mask = clippedmask)
    
    densityOCT17 <- derived(g0.2)
  })
  
  #Double checking where the traps are located  
  
  trapPoints <- traps(oct17)
  par(mfrow = c(1,1), mar = c(1,1,1,1))
  plot(clippedmask, border = 100, ppoly = FALSE)
  gm <- as.data.frame(geom(meio))
  polygon(data.frame(x = gm$x, y = gm$y), col = 'lightgreen', border = NA)
  if (plotSamplingAreas){
    plot(main = "Second Sampling", meio, col = "lightblue")
    plot(trapPoints, detpar = list(pch = 16, cex = 0.8), add = TRUE,
         col = "red")
  }
  
  #Home range results (sigma estimate, lcl exce= lowe95%CI; ucl = upper95%CI):
  # summary(g0.2)$predicted
  
  #Density results (D estimate, lcl = lowe95%CI; ucl = upper95%CI):
  # densityOCT17
  
  # Now we make a table
  cycles <- c(1,2,4,5)
  nCycles <- length(cycles)
  
  blackRatData <- data.table(species = rep("Rattus rattus", times = nCycles), 
                             island = rep("Ilha do Meio", times = nCycles),
                             cycle = cycles,
                             variable = "populationDensity",
                             estimate = c(densityAUG17["D", "estimate"], 
                                          densityOCT17["D", "estimate"],
                                          rep(0, times = 2)),
                             standardError = c(densityAUG17["D", "SE.estimate"], 
                                               densityOCT17["D", "SE.estimate"],
                                               rep(0, times = 2)),
                             LCI = c(densityAUG17["D", "lcl"], 
                                     densityOCT17["D", "lcl"],
                                     rep(0, times = 2)),
                             UCI = c(densityAUG17["D", "ucl"], 
                                     densityOCT17["D", "ucl"],
                                     rep(0, times = 2))
  )
  return(blackRatData)
}