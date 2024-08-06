loadParametersTable <- function(){
  paramTab <- data.table(variable = c("lam(.)", "lam(landscapeCover)",
                                        "gamma(.)", "gamma(totalRainfallThreeMonths)", 
                                      "gamma(TSE)", "gamma(totalRainfallThreeMonths+TSE)",
                                        "p(.)", "p(totalRainfallThreeMonths)", "p(observerID)", 
                                      "p(totalRainfallThreeMonths+observerID)", "p(moonPhase)", "p(totalRainfallThreeMonths+moonPhase)",
                                        "iota(.)", "iota(totalRainfallThreeMonths)", "iota(TSE)", "iota(totalRainfallThreeMonths+TSE)"),
                           value = c(~1, ~landscapeCover,
                                     ~1, ~scale(totalRainfallThreeMonths), ~scale(TSE), ~scale(totalRainfallThreeMonths)+scale(TSE),
                                     ~1, ~scale(totalRainfallThreeMonths), ~observerID, ~scale(totalRainfallThreeMonths)+observerID, ~scale(moonPhase), ~scale(totalRainfallThreeMonths)+scale(moonPhase),
                                     ~1, ~scale(totalRainfallThreeMonths), ~scale(TSE), ~scale(totalRainfallThreeMonths)+scale(TSE)))
  return(paramTab)
}

spParams <- function(islandShort, spShort, addImmigrationFor){
  if (islandShort == "Rata"){
    ########################## RATA
    if (spShort == "elaenia"){
      #### ELAENIA RATA ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(totalRainfallThreeMonths)")
      p <- "p(.)"
      iota <- "iota(.)"
    }
    if (spShort == "mabuia"){
      #### MABUIA RATA ####
      lambda <- c("lam(.)", "lam(landscapeCover)")
      gamma <- "gamma(.)" # no "gamma(totalRainfallThreeMonths)" due to plasticity in feeding!
      p <- c("p(.)", "p(observerID)")
      iota <- "iota(.)"
    }
    if (spShort == "crab"){
      #### CRAB RATA ####
      lambda <- c("lam(.)", "lam(landscapeCover)")
      gamma <- c("gamma(.)")#, "gamma(totalRainfallThreeMonths)")
      p <- c("p(.)")#, "p(totalRainfallThreeMonths)", "p(moonPhase)","p(totalRainfallThreeMonths+moonPhase)")
      iota <- "iota(.)"
    }
    if (spShort == "maskedBooby"){
      #### BOOBIES RATA ####
      lambda <- "lam(.)"
      gamma <- "gamma(.)"
      p <- "p(.)"
      iota <- c("iota(.)")
    }
  } else {
    ########################## MEIO
    if (spShort == "elaenia"){
      #### ELAENIA MEIO ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(totalRainfallThreeMonths)", "gamma(TSE)", "gamma(totalRainfallThreeMonths+TSE)")
      p <- "p(.)"
      iota <- c("iota(.)")
    }
    if (spShort == "mabuia"){
      #### MABUIA MEIO ####
      lambda <- c("lam(.)", "lam(landscapeCover)")
      gamma <- c("gamma(.)", "gamma(TSE)")
      p <- c("p(.)", "p(observerID)")
      iota <- "iota(.)"
    }
    if (spShort == "crab"){
      #### CRAB MEIO ####
      lambda <- c("lam(.)", "lam(landscapeCover)")
      gamma <- c("gamma(.)"), "gamma(totalRainfallThreeMonths)", "gamma(TSE)", "gamma(totalRainfallThreeMonths+TSE)")
      p <- c("p(.)"), "p(totalRainfallThreeMonths)", "p(moonPhase)","p(totalRainfallThreeMonths+moonPhase)")
      iota <- "iota(.)"
    }
    if (spShort == "maskedBooby"){
      #### BOOBIES MEIO ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(TSE)")
      p <- "p(.)"
      iota <- c("iota(.)")
      if (spShort %in% addImmigrationFor){
        iota <- c(iota, "iota(TSE)") 
      }
    }
  }
  return(list(lambda = lambda,
              gamma = gamma,
              p = p,
              iota = iota))
}

