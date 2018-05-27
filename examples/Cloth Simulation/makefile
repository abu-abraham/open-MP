default: kernel_main_basic

kernel_main_basic:
	module load papi;
	gcc -c  kernel_main_basic.c -o cloth_code.o -pg -lm 
	gcc kernel_main.c -o kernel_main_basic cloth_code.o -pg -lm -l papi   


kernel_main_opt:
	module load papi;
	gcc -c  kernel_main_opt.c -o cloth_code.o -pg -lm 
	gcc kernel_main.c -o kernel_main_opt cloth_code.o -pg -lm -l papi  

kernel_main_sse:
	module load papi;
	gcc -c  kernel_main_sse.c -o cloth_code.o -pg -lm 
	gcc kernel_main.c -o kernel_main_sse cloth_code.o -pg -lm -l papi  

kernel_main_opt_block:
	module load papi;
	gcc -c  kernel_main_omp_block.c -o cloth_code.o -pg -lm -fopenmp
	gcc kernel_main.c -o kernel_main_omp cloth_code.o -pg -lm -l papi  -fopenmp

kernel_main_omp:
	module load papi;
	gcc -c  kernel_main_omp.c -o cloth_code.o -pg -lm -fopenmp
	gcc kernel_main.c -o kernel_main_omp cloth_code.o -pg -lm -l papi  -fopenmp

