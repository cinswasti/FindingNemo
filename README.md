# FindingNemo
Compilation of R functions to produce ultra-long sequencing metrics and graphs from ONT platforms, but can be used for any nanopore sequencing output (including ligation based).   
After copying/downloading into a folder, run `source(functionName)` to read and execute the function.  
Required arguments can be called by `args(functionName)` command.   
Annotations of the arguments are given inside every function file, below are some examples.  

## compileOccupancy
Calculate relative occupancies of the three types of available pores (i.e. adapter, pore, and strand)
```r
compileOccupancy(dutyTimeFile = "path/to/duty_time_file/or/pore_activity_file", # different chemistry version calls the file differently
                 workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
                 saveFile = "test_occupancy",
                 startPoint = 11, stopPoint = NULL, runStatus = "run01")
```
Output `test_occupancy.csv` file is saved in `workDir`.

## lengthHistogram
Create a histogram of read length distribution and N50
```r
lengthHistogram(seqSummary = "path/to/sequencing_summary_file",
                workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
                runID = "run01", endTime = 60, yInt = 10, breaksStart = 2, breaksVal = 1)
```
Output `run01_timed_histogram.png` is saved in the same folder with the `sequencing_summary_file.txt`.

## metricsTable
Calculate sequencing metrics, including normalised yield (i.e. yield per pore over a period of data collection time)
```r
metricsTable(seqSummary = "path/to/sequencing_summary_file",
             workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
             runStatus = "run01", stopPoint = 1440,
             saveDir = "path/to/output/directory",
             saveFile = "test_metrics")
```
Output `test_metrics.csv` is saved in `saveDir`.

## yieldStack
Create a stackplot of yield based on read length categories
```r
yieldStack(seqSummary = "path/to/sequencing_summary_file",
           workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
           runID = "run01", platform = "P24", endTime = 120, ystack = 5)
```
Output `run01_timed_stackplot.png` is saved in the same folder with the `sequencing_summary_file.txt`.  

