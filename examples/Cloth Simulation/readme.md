# ReadMe

## To Run 

A makefile is provided which could be executed by typing make followed by the operation. By default, the initial step of ensuring only correct componenets of velocity are set to 0 is executed.

All files are compiled by loading PAPI modules and executed with *gprof*

#### Run Step 1 (Initial Step ) 

> make kernel_main_basic
> ./kernel_main_basic

#### Run Step 2 (Modified Code) 

> make kernel_main_opt
> ./kernel_main_opt

#### Run Step 3 (SSE Vectorization - Only ) 

> make kernel_main_sse
> ./kernel_main_sse

#### Run Step 4 (OpenMP +SSE ) 

> make kernel_main_omp
> ./kernel_main_omp

#### Run Step 5 (OpenMP + SSE + block-wise execution ) 

> make kernel_main_omp_block
> ./kernel_main_omp_block


## To see profiling results

To see profiling results, just execute
> gprof kernel_main gmon.out

*only the last run values are displayed*

