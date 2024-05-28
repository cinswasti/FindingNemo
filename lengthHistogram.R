### Last edited 24/05/2024
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
library(tools)
#
### 1 - Preparing data table for histogram
lengthHistogram <- function(seqSummary = NULL, runID, maxRead = 2000000, xInt = 1000, yInt = 750, 
                            workDir, startTime = 0, endTime = NULL, breaksStart = 25, breaksVal = 25, binSize = 1e+4)
{
  ### Arguments 
  ### seqSummary: path to sequencing_summary file
  ### runID     : runID for graph title 
  ### maxRead   : an integer for maximum read-length to bin all reads (in bp); don't put value < 1e+6 as it will change binning visualization
  ### xInt, yInt: max limit of x- and y- axis, in kb and Mb respectively
  ### workDir   : path to working directory; preferably the parent folder where all the runs are
  ### startTime : minute to start data inclusion
  ### endTime   : mintue to end data inclusion
  ### breaksStart : first y-axis value shown
  ### breaksVal : y-axis interval 
  ### binSize   : histogram bin size (bp)
  #
  setwd(workDir)
  message("------- Preparing Table -------")
  tmp <- read.table(file = seqSummary, header = TRUE, sep = "\t") %>% 
      select(sequence_length_template, start_time, duration) %>% 
      mutate(minute = floor((start_time + duration)/60)) %>% 
      arrange(minute) %>% 
      filter(minute >= startTime & minute <= endTime) %>% 
      rename(bp_length = 1)
    lengthTable <- tmp %>% 
      select(bp_length) %>% 
      arrange(bp_length) %>% 
      group_by(bp_length) %>% 
      mutate(bin = cut(bp_length, breaks = seq(0, maxRead, by = binSize), include.lowest = FALSE, 
                       dig.lab = 0, labels = seq(10, maxRead/1000, by = 10), ordered_result = TRUE),
             bp_yield = sum(bp_length)) %>% 
      unique() %>% 
      group_by(bin) %>% 
      mutate(bin_yield = sum(bp_yield)) %>% 
      select(bin, bin_yield) %>% 
      unique()
  message("------- Table of lengths is created -------")
  #
  N50 <- tmp %>%
      filter(minute <= endTime) %>%
      arrange(desc(bp_length)) %>%
      mutate(bp_sum = sum(as.numeric(bp_length)), bp_cumsum = cumsum(as.numeric(bp_length))) %>%
      filter(bp_cumsum >= bp_sum * 0.5) %>%
      filter(bp_cumsum == min(bp_cumsum)) %>% 
      select(bp_length)
  message("------- N50 of minute ", paste0(endTime), " is ", paste0(N50/1000), " kb -------")
  #
### 2 - Building histogram of read-length bins vs Yield
theme_set(theme_minimal())
cols <- c("1" = "#009E73", "2" = "#CC79A7", "3" = "#0072B2", "4" = "#000000", "5" = "#999999")
histogramLength <- lengthTable %>%
  mutate(bin = as.numeric(as.character(bin))) %>%
  mutate(group = ifelse(bin <= 100, "1",
                        ifelse(bin > 100 & bin <= 200, "2",
                               ifelse(bin > 200 & bin <= 500, "3",
                                      ifelse(bin > 500 & bin < 1000, "4",
                                             ifelse(bin >= 1000, "5", NA)))))) %>%
  ggplot() +
  geom_col(aes(x = bin, y = bin_yield/1e+9, fill = group)) +
  labs(x="Read Length (kb)", y="Yield (Gb)", title = runID) +
  coord_cartesian(xlim = c(0,xInt), ylim = c(0,yInt)) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20), axis.title = element_text(size = 24),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.text = element_text(size = 18, margin = margin(t = 6, r = 0, b = 6, l = 0)),
        legend.title = element_text(size = 20), legend.spacing.y = unit(20, "points"),
        legend.box.spacing = unit(28, "points"),
        title = element_text(size = 24)) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(from = breaksStart, to = yInt, by = breaksVal))) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(from = 250, to = maxRead, by = 250))) +
  geom_vline(data = N50, aes(xintercept = bp_length/1000), color = "black", linewidth = 1.5, linetype = "dashed") +
  scale_fill_manual(values = cols, name = "Read Length",
                    labels = c("Mb <100kb", "Mb 100-199kb", "Mb 200-499kb", "Mb 500-999kb", "Mb >=1000kb"))
# 
## 3 - Outputs
path = dirname(seqSummary)
setwd(path)
plot(histogramLength)
ggsave(paste0(runID, "_timed_histogram.png"), units = "mm", width = 297, height = 210, bg = "white")
setwd(workDir)
}
### END of function