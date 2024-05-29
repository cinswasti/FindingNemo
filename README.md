# FindingNemo
Compilation of R functions to produce ultra-long sequencing metrics and graphs from ONT platforms  
After copying/downloading into a folder, run `source(functionName)` to read and execute the function.  
Required arguments can be called by `args(functionName)` command.   
Annotations of the arguments are given inside every function file, below are some examples.  

## compileOccupancy
Calculate relative occupancies of three types of available pores (i.e. adapter, pore, and strand)
```r
compileOccupancy(dutyTimeFile = "path/to/duty_time_file/or/pore_activity_file", # different chemistry version calls the file differently
                 workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
                 saveFile = "test_occupancy",
                 startPoint = 11, stopPoint = NULL, runStatus = "run01")
```
Output `test_occupancy.csv` file is saved in `workDir`

## lengthHistogram

```r
lengthHistogram(seqSummary = "path/to/sequencing_summary_file",
                workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
                runID = "run01", endTime = 60, yInt = 10, breaksStart = 2, breaksVal = 1)
```
Output `run01_timed_histogram.png` is saved in the same folder with the `sequencing_summary_file.txt`

## metricsTable
```r
metricsTable(seqSummary = "path/to/sequencing_summary_file",
             workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
             runStatus = "run01",
             saveDir = "path/to/output/directory",
             saveFile = "test_metrics")
```

## yieldStack
```r
```

