#include "myGauss.h"
#include <math.h>
#include <omp.h>
#include <time.h>

#define PI 3.141592653589793238462

myGauss::myGauss(double paramValueVal){
    paramValue = paramValueVal;
    oneOverTwoSigSq = 1.0/(2*paramValue*paramValue);
    sqrtTwoPi = sqrt(2*PI);
};

myGauss::myGauss(){
    paramValue = 1.0;
    oneOverTwoSigSq = 1.0/(2*paramValue*paramValue);
    sqrtTwoPi = sqrt(2*PI);
};

myGauss::~myGauss(){
};

double myGauss::getParameter(){
    return paramValue;
};

void myGauss::setParameter(double newSigma){
    paramValue = newSigma;
    oneOverTwoSigSq = 1.0/(2*paramValue*paramValue);
};

double myGauss::evaluate(double value){
    return exp(-value*value*oneOverTwoSigSq);
};

double myGauss::evaluate(double * dataSet, int dataLength){
    double NLL = 0;
    double norm = 1.0/normValue();
    #pragma omp parallel for default(none) shared(dataSet,dataLength,norm) reduction(+:NLL)
    for(int i=0; i<dataLength; i++){
       double evaluated = exp(-dataSet[i]*dataSet[i]*oneOverTwoSigSq);
       evaluated*=norm;
       NLL+=log(evaluated);
    }
    return -NLL;
};

double myGauss::normValue(){
    return paramValue*sqrtTwoPi;
};

double * myGauss::generateData(int length){
    double * output = new double[length];
    srand(time(NULL));
    for(int i=0; i<length; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;
        //Now doing the box muller method
        output[i] = sqrt(-2.0*log(rand1))*cos(2*PI*rand2)*paramValue;
    }
    return output;
};
