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
metricsTable <- function(seqSummary, workDir, runStatus, startPoint = 11, stopPoint = NULL, saveDir, saveFile)
  {
  ## seqSummary = path to sequencing_summary file
  ## workDir    = path to working directory; preferably the parent folder where all the runs are
  ## runStatus  = run label or ID
  ## startPoint = minute to start data inclusion, 11 to exclude the first 10 minutes of mux scan
  ## stopPoint  = minute to stop data inclusion
  ## saveDir    = directory to save output file
  ## saveFile   = name of output file, also include extension (CSV)
  #
  setwd(workDir)
  tab <- read.table(file = seqSummary, header = T, sep = "\t") %>%
      select(-1,-2,-4,-19) %>%
      rename(start = start_time, bp_length = sequence_length_template) %>%
      mutate(end = round(start + duration, 2), minute = floor(end/60)) %>%
      filter(minute >= startPoint & minute <= stopPoint) %>%
      group_by(minute) %>% arrange(.by_group = TRUE)
      #
## 1 - longread yield
long <- tab %>% filter(bp_length >= 100000)
compile <- NULL
compile  <- tibble(unique(long$experiment_id), unique(long$sample_id), sum(long$bp_length))
message("Part 1 - Long read yield created \nPart 2 - N50 calculations started")
#
## 2 - N50 calculation
N50 <- tab %>%
  ungroup(minute) %>%
  arrange(desc(bp_length)) %>%
  mutate(bp_sum = sum(as.numeric(bp_length)), bp_cumsum = cumsum(as.numeric(bp_length))) %>%
  filter(bp_cumsum >= bp_sum * 0.5) %>%
  filter(bp_cumsum == min(bp_cumsum)) %>%
  rename(N50 = bp_length) %>%
  select(N50, bp_sum, bp_cumsum, minute)
message(">>> ", runStatus, " --- N50 at minute ", stopPoint, " is ", N50$N50, " bp")
compile <- add_column(compile, N50$N50, runStatus)
#
## 3 - total yield
compile <- add_column(compile, sum(tab$bp_length))
message("Part 3 - Total yield added to table")
#
## 4 - sequencing pore number
pore <- tab %>%
  ungroup(minute) %>%
  select(channel, mux) %>%
  arrange(channel, mux) %>%
  unique() %>%
  summarise(nrow(.))
compile <- add_column(compile, pore)
message("Part 4 - Active pore number calculated")
#
## 5 - passed bases and reads calculations
pass <- tab %>%
  rename(pass = passes_filtering) %>%
  filter(pass == TRUE)
compile <- add_column(compile, sum(pass$bp_length), nrow(pass), nrow(tab)) %>% 
  mutate(run_minutes = stopPoint, 
         ratio_longread = sum(long$bp_length)/sum(tab$bp_length),
         ratio_shortread = 1 - ratio_longread,
         norm_yield = sum(tab$bp_length)/pore) %>% 
  rename("normalised_yield(bp/pore)" = norm_yield)
message("Part 5 - Passed bases/reads calculated and more parameters added")
#
names(compile)[1:10] <- c("runID", "run_info", "UL_yield(bp)", "finalN50(bp)", "runStatus", 
                          "total_yield(bp)", "seq_pores", "pass_yield(bp)", "pass_reads", "total_reads")
setwd(saveDir)
write.table(compile, file = paste0(saveFile, ".csv"), append = TRUE, row.names = FALSE,
            col.names = !file.exists(saveFile), sep = "\t", quote = FALSE)
message("End part - Table is written to saved directory")
setwd(workDir)
}
#
## END of function
