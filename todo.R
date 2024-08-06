# TODO
# [ ] Run Elaenia for both islands
# [ ] Evaluate the results
# [ ] Run posthoc and make plots
# [ ] Write the results / discussion
# [ ] Cleanup code
# [ ] Publish code as release
# [ ] Send to co-authors
# [ ] Submit paper

# TODO Add test to check if data matches weather and wildlife
 
shortSpIsland <- c("mabuiaMeio", "elaeniaMeio", "crabMeio", "maskedBoobyMeio",
                   "mabuiaRata", "elaeniaRata", "crabRata", "maskedBoobyRata")

folderID <- "14eYtGSpG5OZrsc9xUy88YPQehSI0Rw7p"
drive_upload("~/projects/analysisIlhaDoMeio/outputs.7z", as_id(folderID))


# 7z a outputs ./outputs/*