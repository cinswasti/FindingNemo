# FindingNemo
Compilation of R functions to produce ultra-long sequencing metrics and graphs from ONT platforms  
After copying/downloading into a folder, run `source(functionName)` to read and execute the function.  
Required arguments can be called by `args(functionName)` command.   
Annotations of the arguments are given inside each function itself, below are some examples.  
## compileOccupancy
```r
compileOccupancy(dutyTimeFile = "path/to/duty_time_file/or/pore_activity_file", # different chemistry version calls the file differently
                 workDir = "path/to/working/directory", # better to use the location of the sequencing output folder
                 saveFile = "test_occupancy",
                 startPoint = 11, stopPoint = NULL, runStatus = "run01")
```
