# OpenMP

OpenMP is an API for writting shared memory parallel applications in C, C++ and FORTRAN. This repository contains a slight introduction to openMP and some example codes. 

The usage of OpenMP enables application developers to easly tke advantage of multiple CPU's to boost performance of the sequential code. Compared to other parallel execution models such as PThread ( POSIX Threads ), OpenMP resembles serial code and is much easier to implement.

> **Note:** To use OpenMP,  we have to include (import) omp.h

Being a thread centric model, where each thread has a well defined task assignment, openMP allows nesting of threads. Threads in a team has ranks from 0 to N, with 0 being the master thread. In a parallel execution, master thread joins the team and contributes to the parallel computations. 


## Writing parallel code with openMP

OpenMP allows the use of directives which apply to the code block that follows. An example would be, 

>#pragma omp parallel {
  //code block to be executed in parallel
}

We can restrict the number of threads, the visibility of variables used and various other attributes by adding *clauses*. 

Some of the directives used in openMP is mentioned below: 

> #pragma omp master  **->** ensures that only master thread executes the following code

> #pragma omp single **->** ensures that only a single thread executes the following code. If the executing thread needs to be specified, we can use additional clause. 

Data sharing between variables within a block executed in parallel can be done using options :  private, shared, default, firstprivate etc... The following table illustrates the visibility of each option. 

```
| 			    | Parallel | Sections |  Task  |  Single | for
| :-----------: |:--------:|:--------:|:------:|:-------:|:----:
| private       |    [x]   |   [x]    |   [x]  |    [x]  | [x]
| shared        |    [x]   |          |   [x]  |         | 
| default       |    [x]   |          |   [x]  |         | 
| firstprivate  |    [x]   |   [x]    |   [x]  |    [x]  | [x]
| lastprivate   |          |   [x]    |        |         | [x]
```

**private** ensures that the variables decalred are private for to each task. Whereas **shared** defines variables to be sahed by all tasks.  **firstprivate**  is used when we need to make variable private and intialize it with value of item when parallel block is encountered. It could be viewed as a pass-by-value mechanism. Opposite to firstprivate, **lastprivate** is used to transfers values from the parallel region to the outside context. 

An example snippet: 

   ``` 
   omp_set_num_threads(3); //execute 3 threads 
#pragma omp parallel private(i)
    {
        printf("value of i is %d \n", i);
        i = i+1;
    }
```

Variables declared inside the main and outside the code block is shared, whereas, the ones initialized inside is local. The above mentioned table helps in understanding how variables declared inside and outside the code block is intepreted. 



### Synchronization directives

A **reduction** clause can be added to the parallel directive. This specifies that the final values of certain variables are combined using the specified operation  at the end of the parallel region.  

   ``` 
  #pragma omp parallel default(shared) reduction(+:i) reduction(*:j) 
         {
            i = omp_get_thread_num()+1;
            j = omp_get_thread_num()+1;
            printf("Thread %i: i=%i,j=%i\n",omp_get_thread_num(),i,j,k);
         }

         printf("After parallel region: i=%i,j=%i\n",i,j);
```

**barrier** directive enables us to synchronize threads. It ensures that all threads are complete before executing the following block. 
**critical** directive implements a mutaul exclution to the code block that follows. In effect, it serializes the parallel region by ensuring that only one thread runs the specified block. 
**atomic** directive makes sure that no 2 threads update the memory at the same time. Read operations are fine, but no write. 
**flush** is used to flush regester values to memory. 


Example code for each can be found in  [folder](examples/)

