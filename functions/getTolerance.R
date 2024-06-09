getTolerance <- function(sp, isld){
  dt <- data.table::data.table(species = c(rep("maskedBooby", times = 2),
                          rep("crab", times = 2),
                          rep("elaenia", times = 2),
                          rep("mabuia", times = 2)),
                   island = rep(c("Meio", "Rata"), times = 4),
                   tolerance = c(0.1, 0.8,
                                 0.1, 0.8,
                                 0.5, 0.5,
                                 0.1, 0.1))
 return(dt[species == sp & island == isld, tolerance])
}