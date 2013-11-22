#include "myExp.h"
#include <omp.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <cilk/cilk.h>
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

void myExp::generateData(int dataLength, double * p){
    p = new double[dataLength];
    srand(time(NULL));
    int i = 0;
    while(i<dataLength){ 
        double rand1 = ((double)rand())/RAND_MAX;
            p[i] = log(rand1)/paramValue;
            i++;
        
    }
};

double myExp::evaluate(_Cilk_shared double *  dataSet, int dataSetLength, int threads){
    double NLL = 0;
    double norm = 1.0/normValue();
    omp_set_num_threads(threads);
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
