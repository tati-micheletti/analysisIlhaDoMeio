loadParametersTable <- function(){
  paramTab <- data.table(variable = c("lam(.)", "lam(cover)",
                                        "gamma(.)", "gamma(JulianDate)", "gamma(TSE)", "gamma(JulianDate+TSE)",
                                        "p(.)", "p(JulianDate)", "p(observerID)", "p(JulianDate+observerID)", "p(moonPhase)", "p(JulianDate+moonPhase)",
                                        "iota(.)", "iota(JulianDate)", "iota(TSE)", "iota(JulianDate+TSE)"),
                           value = c(~1, ~cover,
                                     ~1, ~scale(JulianDate), ~scale(TSE), ~scale(JulianDate)+scale(TSE),
                                     ~1, ~scale(JulianDate), ~observerID, ~scale(JulianDate)+observerID, ~moonPhase, ~scale(JulianDate)+moonPhase,
                                     ~1, ~scale(JulianDate), ~scale(TSE), ~scale(JulianDate)+scale(TSE)))
  return(paramTab)
}

spParams <- function(islandShort, spShort){
  if (islandShort == "Rata"){
    ########################## RATA
    if (spShort == "elaenia"){
      #### ELAENIA RATA ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(JulianDate)")
      p <- "p(.)"
      iota <- c("iota(.)", "iota(JulianDate)")
    }
    if (spShort == "mabuia"){
      #### MABUIA RATA ####
      lambda <- c("lam(.)", "lam(cover)")
      gamma <- c("gamma(.)", "gamma(JulianDate)")
      p <- c("p(.)", "p(JulianDate)", "p(observerID)", "p(JulianDate+observerID)")
      iota <- "iota(.)"
    }
    if (spShort == "crab"){
      #### CRAB RATA ####
      lambda <- c("lam(.)", "lam(cover)")
      gamma <- c("gamma(.)", "gamma(JulianDate)")
      p <- c("p(.)", "p(JulianDate)", "p(moonPhase)","p(JulianDate+moonPhase)")
      iota <- "iota(.)"
    }
    if (spShort == "maskedBooby"){
      #### BOOBIES RATA ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(JulianDate)")
      p <- "p(.)"
      iota <- c("iota(.)", "iota(JulianDate)")
    }
  } else {
    ########################## MEIO
    if (spShort == "elaenia"){
      #### ELAENIA MEIO ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(JulianDate)", "gamma(TSE)", "gamma(JulianDate+TSE)")
      p <- "p(.)"
      iota <- c("iota(.)", "iota(JulianDate)", "iota(TSE)", "iota(JulianDate+TSE)")
    }
    if (spShort == "mabuia"){
      #### MABUIA MEIO ####
      lambda <- c("lam(.)", "lam(cover)")
      gamma <- c("gamma(.)", "gamma(JulianDate)", "gamma(TSE)", "gamma(JulianDate+TSE)")
      p <- c("p(.)", "p(JulianDate)", "p(observerID)", "p(JulianDate+observerID)")
      iota <- "iota(.)"
    }
    if (spShort == "crab"){
      #### CRAB MEIO ####
      lambda <- c("lam(.)", "lam(cover)")
      gamma <- c("gamma(.)", "gamma(JulianDate)", "gamma(TSE)", "gamma(JulianDate+TSE)")
      p <- c("p(.)", "p(JulianDate)", "p(moonPhase)","p(JulianDate+moonPhase)")
      iota <- "iota(.)"
    }
    if (spShort == "maskedBooby"){
      #### BOOBIES MEIO ####
      lambda <- "lam(.)"
      gamma <- c("gamma(.)", "gamma(JulianDate)", "gamma(TSE)", "gamma(JulianDate+TSE)")
      p <- "p(.)"
      iota <- c("iota(.)", "iota(JulianDate)", "iota(TSE)", "iota(JulianDate+TSE)")
    }
  }
  return(list(lambda = lambda,
              gamma = gamma,
              p = p,
              iota = iota))
}
