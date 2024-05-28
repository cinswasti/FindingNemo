### Last edited 28/05/2024
# Setting up environment
### If a package is installed, it will be loaded. If any are not, the missing 
### package(s) will be installed from CRAN and then loaded.
### (copied from https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/)
### First specify the packages of interest
packages = c("dplyr", "tidyverse")
#
### Now load or install & load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#
## 1 - Calculate total occupancy of the 3 states (strand, adapter and pore)
compileOccupancy <- function(dutyTimeFile, startPoint = 0, stopPoint = NULL, runStatus, workDir, saveFile)
  {
  ## dutyTimeFile = duty file (or pore activity file in the newer MinKNOW version)
  ## startPoint = n-th minute to start, 11 to exclude the first mux scan
  ## stopPoint  = n-th minute to stop, default NULL to process till the end
  ## workDir    = pathway to working directory where run folders are located (parent directory)
  ## runStatus  = run name or ID
  ## saveFile   = name of output file, file extension included (.csv), output to workDir
  #
  setwd(workDir)
  tab_3states <- read.csv(dutyTimeFile, header = TRUE) %>%
    rename(state = 1, minutes = 2, samples = 3) %>%
    group_by(minutes) %>%
    arrange(.by_group = TRUE) %>%
    mutate(samples_total = sum(samples), samples_portion = samples/samples_total*100) %>%
    #### calculating the total available pores
    filter(state %in% c("adapter", "pore", "strand")) %>%
    mutate(available_samples = sum(samples), available_samples_portion = samples/available_samples*100) %>%
    filter(available_samples_portion != "NaN")
  message("Table is extracted from input file")
#
## 2 - Calculate relative adapter occupancy
adapt <- tab_3states %>%
  filter(state == "adapter") %>%
  group_by(minutes) %>%
  mutate(sum_occupancy = sum(available_samples_portion)) %>%
  select(minutes, sum_occupancy) %>% unique()
if(is.null(stopPoint))
{
  avgOccupancyAdapt <- mean(adapt$sum_occupancy)
  message("Average occupancy of the adapter has been calculated = ", round(avgOccupancyAdapt, digits = 2), "%")
} else {
  adapt <- adapt %>%
    filter(minutes >= startPoint & minutes <= stopPoint) 
  avgOccupancyAdapt <- mean(adapt$sum_occupancy)
  message(runStatus, " --- Average occupancy of adapter until minute ", stopPoint, " has been calculated = ", 
          round(avgOccupancyAdapt, digits = 2), "%")
}
#
## 3 - Calculate relative strand occupancy
strand <- tab_3states %>%
  filter(state == "strand") %>%
  group_by(minutes) %>%
  mutate(sum_occupancy = sum(available_samples_portion)) %>%
  select(minutes, sum_occupancy) %>% unique()
if(is.null(stopPoint))
{
  avgOccupancyStrand <- mean(strand$sum_occupancy)
  message("Average occupancy of the strand has been calculated = ", round(avgOccupancyStrand, digits = 2), "%")
} else {
  strand <- strand %>%
    filter(minutes >= startPoint & minutes <= stopPoint) 
  avgOccupancyStrand <- mean(strand$sum_occupancy)
  message(runStatus, " --- Average occupancy of strand until minute ", stopPoint, " has been calculated = ", 
          round(avgOccupancyStrand, digits = 2), "%")
}
#
## 4 - Calculate relative pore occupancy
pore <- tab_3states %>%
  filter(state == "pore") %>%
  group_by(minutes) %>%
  mutate(sum_occupancy = sum(available_samples_portion)) %>%
  select(minutes, sum_occupancy) %>% unique()
if(is.null(stopPoint))
{
  avgOccupancyPore <- mean(pore$sum_occupancy)
  stopPoint <- max(pore$minutes)
  message("Average occupancy of the pore has been calculated = ", round(avgOccupancyPore, digits = 2), 
          "\n Run time is ", stopPoint, " minutes")
} else {
  pore <- pore %>%
    filter(minutes >= startPoint & minutes <= stopPoint) 
  avgOccupancyPore <- mean(pore$sum_occupancy)
  message(runStatus, " --- Average occupancy of pore until minute ", stopPoint, " has been calculated = ", 
          round(avgOccupancyPore, digits = 2), "%")
}
## 5 - Compile output table
setwd(workDir)
compile <- tibble(runStatus, stopPoint, avgOccupancyAdapt, avgOccupancyStrand, avgOccupancyPore)
names(compile) <- c("runID", "end_minute", "adapter%", "strand%", "pore%")
write.table(compile, file = paste0(saveFile,".csv"), append = TRUE, row.names = FALSE,
                       col.names = !file.exists(saveFile), sep = "\t", quote = FALSE)
}
#
## END of function
