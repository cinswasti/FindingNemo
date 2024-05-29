# FindingNemo
Compilation of R functions to produce ultra-long sequencing metrics and graphs from ONT platforms
After copying/downloading into a folder, run `source(functionName)` to read and execute the function.
Required arguments can be called by `args(functionName)` command.
Annotations of the arguments are given inside each function itself, below are some examples.
## compileOccupancy
```r
compileOccupancy(dutyTimeFile = "path/to/duty_time_file/or/pore_activity_file",
                 workDir = "path/to/directory/where/sequencing/output/folder/is",
                 saveFile = "test_occupancy",
                 startPoint = 11, stopPoint = NULL, runStatus = "run01")
```
