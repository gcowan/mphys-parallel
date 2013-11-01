#include "myGauss.h"
#include <math.h>

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

double myGauss::normValue(){
    return sigma*sqrtTwoPi;
};
