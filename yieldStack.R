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
library(tools)
#
### 1 - Table preparation
yieldStack <- function(seqSummary, runID, platform = NULL, ystack = 10, workDir, 
                       startTime = 0, endTime = NULL, binSize = 1e+4, maxRead = 2e+6)
{
  ### Arguments
  ### seqSummary: sequencing_summary file; seqTable : if there is a table already made
  ### runID     : runID for graph title
  ### platform  : sequencing platform used (MinION, GridION,PromethION, etc.)
  ### ystack    : maximum y-axis value (Gb)
  ### workDir   : path to working directory; preferably the parent folder where all the runs are
  ### startTime : minute to start data inclusion
  ### endTime   : mintue to end data inclusion
  ### binSize   : histogram bin size (bp)
  ### maxRead   : an integer for maximum read-length to bin all reads (in bp); don't put value < 1e+6 as it will change binning visualization
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
### 2 - Building stacked plot graph of read length groups vs total yield
  cols <- c("1" = "#009E73", "2" = "#CC79A7", "3" = "#0072B2", "4" = "#000000", "5" = "#999999")
  stackplotLength <- lengthTable %>%
    mutate(bin = as.numeric(as.character(bin))) %>%
    mutate(group = ifelse(bin <= 100, "1",
                          ifelse(bin > 100 & bin <= 200, "2",
                                 ifelse(bin > 200 & bin <= 500, "3",
                                        ifelse(bin > 500 & bin < 1000, "4",
                                               ifelse(bin >= 1000, "5", NA)))))) %>%
    mutate(category = platform) %>%
    ggplot(aes(x = category, y = (bin_yield/1e+9), fill = group)) +
    geom_col(position = "stack", width = 0.2) +
    labs(x=NULL, y="Yield (Gb)", title = runID) +
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          axis.text = element_text(size = 20), axis.title = element_text(size = 24),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          legend.text = element_text(size = 18, margin = margin(t = 6, r = 0, b = 6, l = 0)),
          legend.title = element_text(size = 20), legend.spacing.y = unit(20, "points"),
          title = element_text(size = 24)) +
    scale_x_discrete(expand = expansion(mult = c(0.15, 0.5))) +
    ylim(0, ystack) +
    scale_fill_manual(values = cols, name = "Read Length",
                    labels = c("Mb <100kb", "Mb 100-199kb", "Mb 200-499kb", "Mb 500-999kb", "Mb >=1000kb"))
#
## 3 - Outputs
  path = dirname(seqSummary)
  setwd(path)
  plot(stackplotLength)
  ggsave(paste0(runID, "_timed_stackplot.png"), units = "mm", width = 210, height = 297, bg = "white")
  setwd(workDir)
}
### END of function