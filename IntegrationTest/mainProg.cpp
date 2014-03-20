#define PI 3.141592653589793238462
#define  FREE  alloc_if(0) free_if(1)
#define  REUSE  alloc_if(0) free_if(0)
#define  ALLOC  alloc_if(1) free_if(0)

#pragma offload_attribute(push,target(mic))
double * params;
double * limits;
#include <omp.h>
#include "fit.h"
#include <stdio.h>
#pragma offload_attribute(pop)

//Minuit interface
#include <mkl_vsl.h>




int main(){
    params = (double *)malloc(sizeof(double)*dimensions);
    for(int i=0; i<dimensions; i++){
        params[i] = 0.783;
    }
    limits = (double *)malloc(sizeof(double)*dimensions*2);
    for(int i=0; i<dimensions*2; i+=2 ){
        limits[i] = -5.0;
        limits[i+1] = 5.0;
    }
    int threads = 240;
    double time = 0;
    double timeIn = 0;
    double result = 0;
     #pragma offload_transfer target(mic:0) in(limits : length(dimensions*2) ALLOC)
    for(int i=0; i<10; i++){
        time = omp_get_wtime(); 
       #pragma offload target(mic:0) inout(result,timeIn) in(params : length(dimensions) ALLOC) nocopy(limits : REUSE)
        {
            timeIn = omp_get_wtime();
            result = integrateVegas(limits,threads,params);
            timeIn = omp_get_wtime()-timeIn;
        }
        time = omp_get_wtime()-time;
        printf("That integration took %f seconds Xeon Phi runtime=%f \n",time,timeIn);
    }
    return 0;
}






