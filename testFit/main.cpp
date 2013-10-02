#include "negativeLog.h"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <iostream>
#include <chrono>

//Method to find max value in an array returning its index
int findMax(std::vector<double>);
int findMin(std::vector<double>);


int main(){
    const double PI  =3.141592653589793238462;
    double sigma = 0.783;
    //Creating the timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
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
    start = std::chrono::system_clock::now();
    for(int i=0; i<sigmaLength; i++){
        double sigmaVal = (((double) i)/sigmaLength)*2+0.001;
        negativeLogs.push_back(test.evaluateDataSetSerial(data,length,sigmaVal));
        
    }
    //finding the min value of loglikelihood 
    int minIndex = findMin(negativeLogs);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end-start; 
    std::cout<<"Serial: The value that fits this gaussian is: "<<(((double) minIndex)/sigmaLength)*2+0.001<<" this took: "<<elapsedSeconds.count()<<" seconds"<<std::endl;
    
    //Clearing the vector
    negativeLogs.clear();

    //Creating array to store values of the negative log likelihood 
    start = std::chrono::system_clock::now();
    for(int i=0; i<sigmaLength; i++){
        double sigmaVal = (((double) i)/sigmaLength)*2+0.001;
        negativeLogs.push_back(test.evaluateDataSetParallel(data,length,sigmaVal));
        
    }
    //finding the min value of loglikelihood 
    minIndex = findMin(negativeLogs);
    end = std::chrono::system_clock::now();
    elapsedSeconds = end-start; 
    std::cout<<"Parallel: The value that fits this gaussian is: "<<(((double) minIndex)/sigmaLength)*2+0.001<<" this took: "<<elapsedSeconds.count()<<" seconds"<<std::endl;
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

