#include "negativeLog.h"
#include <stdio>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <climits>
//Method to find max value in an array returning its index
int findMax(std::vector<double>);
int findMin(std::vector<double>);


int main(){
    const double PI  =3.141592653589793238462;
    //generating the normally distributed doubles using the box mueller method
    int length = 10000000;
    double * data = new double[length];
    srand(time(NULL));
    for(int i =0; i<length; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;

        //Now doing the box muller method
        data[i] = arccos(sqrt(-2.0*log(rand1))*cos(2*PI*rand2));
    }
    
    //Creating array to store values of the negative log likelihood 
    int sigmaLength = 1000;
    std::vector<double> negativeLogs;
    negativeLog test;
    for(int i=0; i<sigmaLength; i++){
        double sigmaVal = (((double) i)/sigmaLength)*2+0.001;
        negativeLogs.push_back(test.evaluateDataSet(data,length,sigmaVal));
        
    }
    //finding the min value of loglikelihood 
    int minIndex = findMin(negativeLogs);
    std::cout<<"The value that fits this gaussian is"<<(((double) minIndex)/sigmaLength)*2+0.001<<std::endl;


}

int findMax(vector<double> data){
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

int findMin(vector<double> data){
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

