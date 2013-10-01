#include "gaussianFunc.h"
#include <cmath>

gaussianFunc::gaussianFunc(double sigmaVal){
    sigma = sigmaVal;
};

double gaussianFunc::evaluate(double x){
    double out = exp(-(x*x)/(2*sigma*sigma));
    return out;
};


double gaussianFunc::integral(){
    double Pi  =3.141592653589793238462;
    double out = sigma*sqrt(2*Pi);
    return out;
};

void gaussianFunc::setSigma(double newSigma){
    sigma = newSigma;
};

double gaussianFunc::getSigma(){
    return sigma;
}
