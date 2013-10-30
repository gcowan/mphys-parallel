#include "negativeLog.h"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <stdio.h>
#include <omp.h>

//Method to find max value in an array returning its index
int findMax(std::vector<double>);
int findMin(std::vector<double>);


int main(){
    const double PI  =3.141592653589793238462;
    double sigma = 0.783;
    FILE * runTimes = fopen("runTimes.dat","w");
     // FILE * speedup = fopen("speedup.dat","w");
    //generating the normally distributed doubles using the box mueller method
    int length = 10000000;
    double * data = new double[length];
    srand(time(NULL));
    for(int i =0; i<length; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;

        //Now doing the box muller method
        data[i] = sqrt(-2.0*log(rand1))*cos(2*PI*rand2)*sigma;
    }
    
    //Creating array to store values of the negative log likelihood 
    int sigmaLength = 1000;
    std::vector<double> negativeLogs;
    negativeLog test;
    /*
    double time = omp_get_wtime();
    for(int i=0; i<sigmaLength; i++){
        double sigmaVal = (((double) i)/sigmaLength)*2+0.001;
        negativeLogs.push_back(test.evaluateDataSetSerial(data,length,sigmaVal));
        
    }
    //finding the min value of loglikelihood 
    int minIndex = findMin(negativeLogs);
    time = omp_get_wtime()-time;
    double serialTime = time;
    fprintf(runTimes,"%d %f\n",1,serialTime);
     // fprintf(speedup,"%d %f\n",1,1.0);
    //Clearing the vector
    negativeLogs.clear();
   */ 


    //Doing for two to 500 threads to measure the speedup
    for(int threads=100;threads<500; threads+=2){
	omp_set_num_threads(threads);
    	double time = omp_get_wtime();
        for(int i=0; i<sigmaLength; i++){
            double sigmaVal = (((double) i)/sigmaLength)*2+0.001;
            negativeLogs.push_back(test.evaluateDataSetParallel(data,length,sigmaVal));
        }
        //finding the min value of loglikelihood 
        int minIndex = findMin(negativeLogs);
        time = omp_get_wtime()-time;
	    fprintf(runTimes,"%d %f\n",threads,time);
	 // fprintf(speedup,"%d %f\n",threads,serialTime/time);
        negativeLogs.clear();
    }
    fclose(runTimes);
     // fclose(speedup); 
    return 0;




}

int findMax(std::vector<double> data){
    double maxValue = INT_MIN;
    int maxIndex = INT_MIN;
    for(int i=0; i<data.size(); i++){
        if(data[i]>maxValue){
            maxValue = data[i];
            maxIndex = i;
        }
    }
    return maxIndex;
};

int findMin(std::vector<double> data){
    double minValue = INT_MAX;
    int minIndex = INT_MIN;
    for(int i=0; i<data.size(); i++){
        if(data[i]<minValue){
            minValue = data[i];
            minIndex = i;
        }
    }
    return minIndex;
}

