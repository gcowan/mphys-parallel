#include "myGauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <limits.h>
#include <float.h>

int main(){
    const double PI  =3.141592653589793238462;
    double sigma = 0.783;
    int scanNumber = 1000;
    //generating the normally distributed doubles using the box mueller method
    int length = 10000000;
    double * data = new double[length];
    srand(time(NULL));
    for(int i=0; i<length; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;
        //Now doing the box muller method
        data[i] = sqrt(-2.0*log(rand1))*cos(2*PI*rand2)*sigma;
    }
    double * NLLS = new double[scanNumber];
    FILE * output = fopen("unoffloaded.dat","w");
    myGauss gaussp[SIZE];
    for(int threads = 1; threads<=24; threads+=1){
        double time = omp_get_wtime();
         // #pragma offload target(mic:0) in(data:length(length) ) inout(NLLS:length(scanNumber))
        {
            omp_set_num_threads(threads);
            for(int i=0 ; i<scanNumber; i++){
                double sigmaGuess = i*0.001+0.001;
                gaussp[0].setSigma(sigmaGuess);
                //Calculating the negative loglikeliHood
                NLLS[i]=gaussp[0].evaluate(data,length);
            }
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
    double minSig = minIndex*0.001+0.001;
    time = omp_get_wtime()-time;
    if(minSig!=sigma)
        printf("Somthing has gone wrong %f",minSig);
    else
        printf("Sigma val %f\n",minSig);
        fprintf(output,"%d %f\n",threads,time);

    }
    fclose(output);
    return 0;
}
