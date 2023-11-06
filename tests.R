# Further Tests
toRemove <- c("Island", "Year", "Month", "Day", "Species",
"Radius", "Original_Landscape", "originDate",
"pointID", "transectID")

mabuiaMeio <- DT[Island == "Ilha do Meio" & Species == "Trachylepis atlantica",]
mabuiaMeio <- mabuiaMeio[, (toRemove) := NULL]
y1 <- dcast(mabuiaMeio, site ~ Date, value.var = "Counts")
cMeio <- DT[Island == "Ilha do Meio" & Species == "Johngarthia lagostoma",]
cMeio <- cMeio[, (toRemove) := NULL]
y2 <- dcast(cMeio, site ~ Date, value.var = "Counts")
ccMeio <- DT[Island == "Ilha do Meio" & Species == "Elaenia ridleyana",]
ccMeio <- ccMeio[, (toRemove) := NULL]
y3 <- dcast(ccMeio, site ~ Date, value.var = "Counts")
# ======================================================
mabuiaRata <- DT[Island == "Ilha Rata" & Species == "Trachylepis atlantica",]
mabuiaRata <- mabuiaRata[, (toRemove) := NULL]
y4 <- dcast(mabuiaRata, site ~ Date, value.var = "Counts")
cRata <- DT[Island == "Ilha Rata" & Species == "Johngarthia lagostoma",]
cRata <- cRata[, (toRemove) := NULL]
y5 <- dcast(cRata, site ~ Date, value.var = "Counts")
ccRata <- DT[Island == "Ilha Rata" & Species == "Elaenia ridleyana",]
ccRata <- ccRata[, (toRemove) := NULL]
y6 <- dcast(ccRata, site ~ Date, value.var = "Counts")
y1
y2
y3
y4
y5
y6
