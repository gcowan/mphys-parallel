#include "myGauss.h"
#include <math.h>
#include <omp.h>

myGauss::myGauss(double sigmaVal){
    sigma = sigmaVal;
    oneOverTwoSigSq = 1.0/(2*sigma*sigma);
    sqrtTwoPi = sqrt(2*3.14159265359);
};

myGauss::myGauss(){
    sigma = 1.0;
    oneOverTwoSigSq = 1.0/(2*sigma*sigma);
    sqrtTwoPi = sqrt(2*3.14159265359);
};

myGauss::~myGauss(){
};

double myGauss::getSigma(){
    return sigma;
};

void myGauss::setSigma(double newSigma){
    sigma = newSigma;
    oneOverTwoSigSq = 1.0/(2*sigma*sigma);
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
}

double myGauss::normValue(){
    return sigma*sqrtTwoPi;
};
