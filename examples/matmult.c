#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

typedef double TYPE;
#define MAX_VAL 1000
#define MIN_VAL 1
#define N 400

TYPE** generateRandomMatrix(int dimension, int const_value){

	TYPE** matrix = malloc(dimension * sizeof(TYPE*));
    int i,j;
	for(i=0; i<dimension; i++){
		matrix[i] = malloc(dimension * sizeof(TYPE));
	}

	for(i=0; i<dimension; i++){
		for(j=0; j<dimension; j++){
			matrix[i][j] = (rand() % MAX_VAL + MIN_VAL)*const_value;
		}
	}

	return matrix;
}

void normalMatrixMultiplication(TYPE** A, TYPE** B, TYPE** C, int dimension){

    int i,j,k;
	for(i=0; i<dimension; i++){
		for(j=0; j<dimension; j++){
			for(k=0; k<dimension; k++){
				C[i][j] = C[i][j] + (A[i][k] * B[k][j]);
			}
		}
	}

    //To verify whether the multiplication is computed correctly
    double sum = 0;
    for(i=0; i<dimension; i++){
		for(j=0; j<dimension; j++){
            sum+=C[i][j];
        }
    }
    printf("%lf\n",sum);
}

void parellelMultiplicationSubLoop(TYPE** A, TYPE** B, TYPE** C, int dimension, int i, int j){
    int k;
    for(k=0; k<dimension; k++){
		C[i][j] = C[i][j] + (A[i][k] * B[k][j]);
	}
}

void parallelMatrixMultiplication(TYPE** A, TYPE** B, TYPE** C, int dimension){

    int i,j,k;
    omp_set_num_threads(10);
    #pragma omp parallel for 
	for(i=0; i<dimension; i++){
        #pragma omp parallel for 
		for(j=0; j<dimension; j++){
			parellelMultiplicationSubLoop(A,B,C,dimension,i,j);
		}
	}

     
    double sum =0;
    //To verify whether the multiplication is computed correctly
    #pragma omp barrier
	for(i=0; i<dimension; i++){
		for(j=0; j<dimension; j++){
            sum+=C[i][j];
        }
    }
    printf("%lf\n",sum);
}

void main(){
    TYPE ** A = generateRandomMatrix(N,1);
    TYPE ** B = generateRandomMatrix(N,1);
    TYPE ** C = generateRandomMatrix(N,0);
    normalMatrixMultiplication(A,B,C,N);
    C = generateRandomMatrix(N,0);
    parallelMatrixMultiplication(A,B,C,N);
}