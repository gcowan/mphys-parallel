#include "myFunc.h"
#include "myGauss.h"
#include "myExp.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <cilk/cilk.h>
#pragma offload_attribute(push,_Cilk_shared)
#include <offload.h>
#pragma offload_attribute(pop)

using namespace std;

_Cilk_shared myFunc * _Cilk_shared  funcp;
_Cilk_shared const int length = 1000000;
_Cilk_shared  double  *  _Cilk_shared  data;  
_Cilk_shared int threads;
_Cilk_shared myGauss gauss;
_Cilk_shared double * _Cilk_shared NLLS;

double _Cilk_shared  fit(_Cilk_shared myGauss&  functionp, _Cilk_shared double * NLLS, int scanNumber,  int threads,  _Cilk_shared double *  data,  int length){
    for(int i=0 ; i<scanNumber; i++){
        double sigmaGuess = i*0.001+0.001;
        functionp.setParameter(sigmaGuess);
        //Calculating the negative loglikeliHood
        NLLS[i] = functionp.evaluate(data,length,threads);
    }
    //Finding the minimum 
    int minIndex = INT_MAX;
    double minValue = DBL_MAX;
    for(int i=0; i<scanNumber; i++){
        if(NLLS[i]<minValue){
           minValue = NLLS[i];
            minIndex = i;
        }
    }
    double minSig = minIndex*0.001+0.001;
    return minSig;
};


int main(){
    double sigma = 0.783;
    int scanNumber = 1000;
    funcp = &gauss;
    data =  (_Cilk_shared double * ) _Offload_shared_malloc(sizeof(double)*length);
    gauss.setParameter(sigma);
    gauss.generateData(length,data);
    printf("Finished generating data\n");
    NLLS = (_Cilk_shared double *  ) _Offload_shared_malloc(sizeof(double)*scanNumber);
    FILE * output = fopen("offloaded.dat","w");
    for(threads=24; threads<=244; threads+=10){
        double time = omp_get_wtime();
        double minSig =  _Cilk_offload_to(0) fit(gauss,NLLS,scanNumber,threads,data,length);
        time = omp_get_wtime()-time;
        printf("Sigma val %f\n",minSig);
        fprintf(output,"%d %f\n",threads,time);
    }
    fclose(output);
    return 0;
}
