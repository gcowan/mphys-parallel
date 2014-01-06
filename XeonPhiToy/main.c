#include<omp.h>
#include<stdio.h>
#include<math.h>
#include <limits.h>
#include <time.h>
#include <float.h>
#include <stdlib.h>

int main(){
    //generating gaussian stuff using box muller
    const double PI  =3.141592653589793238462;
    double sigma = 0.783;
    int scanNumber = 1000;
    //generating the normally distributed doubles using the box mueller method
    int length = 10000000;
    double * data = malloc(sizeof(double)*length);
    srand(time(NULL));
    for(int i=0; i<length; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;

        //Now doing the box muller method
        data[i] = sqrt(-2.0*log(rand1))*cos(2*PI*rand2)*sigma;
    }
    double * NLLS = malloc(sizeof(double)*scanNumber);
    FILE * output = fopen("unoffloaded.dat","w");
    for(int threads=1; threads<=24; threads+=1){
    double time = omp_get_wtime();
    //Offloading to the Xeon Phi
  //  #pragma offload target(mic:0) in(data:length(length)) inout(NLLS:length(scanNumber))
    {
        //Iterating over the values of sigma
        omp_set_num_threads(threads);
        for(int i=0; i<scanNumber; i++){
            double sigmaGuess = i*0.001+0.001;
            //Calculating the negative loglikeliHood
            double NLL = 0;
            double oneOverTwoSigSq = 1.0/(2*sigmaGuess*sigmaGuess);
            #pragma omp parallel for default(none) shared(sigmaGuess,length,data,oneOverTwoSigSq) reduction(+:NLL) 
            for(int j=0; j<length; j++){
                double value = exp(-(data[j]*data[j]*oneOverTwoSigSq));
                double normalisation = sigmaGuess*sqrt(2*PI);
                NLL = NLL + log(value/normalisation);
                
            }
            NLLS[i]=-NLL;
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
    if(minSig!=0.783)
        printf("Somthing has gone wrong %f",minSig);
    else
        fprintf(output,"%d %f\n",threads,time);
    
    }
    fclose(output);
    return 0;
}
