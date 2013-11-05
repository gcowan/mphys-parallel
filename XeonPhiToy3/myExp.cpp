#include "myExp.h"
#include <omp.h>
#include <math.h>
#include <time.h>
#include <limits.h>

myExp::myExp(){
    paramValue = -1;
};

myExp::myExp(double newExponent){
    paramValue = newExponent;
};

myExp::~myExp(){
};

void myExp::setParameter(double newExponent){
    paramValue = newExponent;
};

double myExp::getParameter(){
    return paramValue;
};

double * myExp::generateData(int dataLength){
    double * output = new double[dataLength];
    srand(time(NULL));
    int i = 0;
    while(i<dataLength){ 
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;
        if(rand1<=exp(-paramValue*rand2)){
            output[i] = rand1;
            i++;
        }
    }
    return output;
};

double myExp::evaluate(double * dataSet, int dataSetLength){
    double NLL = 0;
    double norm = 1.0/normValue();
    #pragma omp parallel for default(none) shared(norm, dataSet, dataSetLength) reduction(+:NLL)
    for(int i=0; i<dataSetLength; i++){
        double value = exp(paramValue*dataSet[i]);
        value *= norm;
        NLL += log(value);
    };
    return -NLL;
};

double myExp::normValue(){
    return -1.0/paramValue;
};
