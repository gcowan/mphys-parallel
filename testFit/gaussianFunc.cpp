#include "gaussianFunction.h"
#include <cmath>

gaussianFunc::gaussianFunc(double sigmaVal){
    sigma = sigmaVal;
};

double gaussianFunc::evaluate(double x){
    double out = exp(-(x*x)/(2*sigma*sigma));
    return out;
};


double gaussianFunc::integral(){
    double out = sigma*sqrt(2*PI);
    return out;
};

void gaussianFunc::setSigma(double newSigma){
    sigma = newSigma;
};

double gaussianFunc::getSigma(){
    return sigma;
}
