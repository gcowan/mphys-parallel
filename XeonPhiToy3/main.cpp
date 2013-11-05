#include "myFunc.h"
#include "myGauss.h"
#include "myExp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <cilk/cilk.h>

int main(){
    double sigma = -0.783;
    int scanNumber = 1000;
    int length = 10000000;
    int threads = 24;
    myFunc * funcp = new myExp(sigma);
    double * data = funcp->generateData(length);
    FILE * generated = fopen("data.dat","w");
    for(int i=0; i<length; i++){
        fprintf(generated,"%f\n",data[i]);
    }
    fclose(generated);
    double * NLLS = new double[scanNumber];
    FILE * output = fopen("offloaded.dat","w");
    double time = omp_get_wtime();
    omp_set_num_threads(threads);
    for(int i=0 ; i<scanNumber; i++){
        double sigmaGuess = -i*0.001-0.001;
       //_Cilk_offload  
           funcp->setParameter(sigmaGuess);
        //Calculating the negative loglikeliHood
        NLLS[i]=funcp->evaluate(data,length);
    }

    
    //Finding the minimum on the host
    int minIndex = INT_MAX;
    double minValue = DBL_MAX;
    for(int i=0; i<scanNumber; i++){
        if(NLLS[i]<minValue){
            minValue = NLLS[i];
            minIndex = i;
        }
    }
    double minSig = minIndex*-0.001-0.001;
    time = omp_get_wtime()-time;
    if(minSig!=sigma)
        printf("Somthing has gone wrong %f\n",minSig);
    else
        printf("Sigma val %f\n",minSig);
        fprintf(output,"%d %f\n",threads,time);
    fclose(output);
    return 0;
}
