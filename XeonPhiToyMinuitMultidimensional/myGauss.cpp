#include "myGauss.h"
#include <math.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#define PI 3.141592653589793238462
//Generate Gaussian with set of parameters
myGauss::myGauss(int numSigmas,  double *  paramValues){
    paramValue = paramValues;
    dimensions = numSigmas;
    oneOverTwoSigsSq = (_Cilk_shared double *)_Offload_shared_aligned_malloc(sizeof(int)*numSigmas,64);
    for(int i=0; i<numSigmas; i++){
        oneOverTwoSigsSq[i] = 1.0/2*(paramValue[i]*paramValue[i]);
    }
    sqrtTwoPi = sqrt(2*PI);
};
//Generate gaussian with all parameters set to one
myGauss::myGauss(int numSigmas){
    paramValue = (_Cilk_shared double *)_Offload_shared_aligned_malloc(sizeof(int)*numSigmas,64);
    oneOverTwoSigsSq = (_Cilk_shared double *)_Offload_shared_aligned_malloc(sizeof(int)*numSigmas,64);
    dimensions = numSigmas;
    for(int i=0; i<numSigmas; i++){
        paramValue[i]=1.0;
        oneOverTwoSigsSq[i] = 1.0/2*(paramValue[i]*paramValue[i]);
    };
    sqrtTwoPi = sqrt(2*PI);
};


//Deconstructor
myGauss::~myGauss(){
    _Offload_shared_aligned_free(paramValue);
    _Offload_shared_aligned_free(oneOverTwoSigsSq);
};

//Method to get a particular parameter
double myGauss::getParameter(int index){
    return paramValue[index];
};

//method to set a particular parameter
void myGauss::setParameter(int index, double newSigma){
    paramValue[index] = newSigma;
    oneOverTwoSigsSq[index] = 1.0/(2*paramValue[index]*paramValue[index]);
};

//Method to evaluate function for one set of x values
double myGauss::evaluate(double * value){
    double output = exp(-value[0]*value[0]*oneOverTwoSigsSq[0]);
    for(int i=1; i<dimensions; i++){
        output *= exp(-value[i]*value[i]*oneOverTwoSigsSq[i]);
    }
    return output;
};

//Method to evaluate entire dataset using a set number of threads
double myGauss::evaluate( double * dataSet, int dataLength, int threads ){
    double NLL = 0;
    double norm = 1.0/normValue();
    omp_set_num_threads(threads);
    /*int off = _Offload_get_device_number();
    if(off==0)
        printf("On mic\n");
    if(off==-1)
        printf("On host\n");
   */ 
    #pragma omp parallel for default(none) shared(dataSet,dataLength,norm) reduction(+:NLL)
    for(int i=0; i<dataLength*dimensions; i+=dimensions){
       double evaluated = evaluate(dataSet+i);
       evaluated*=norm;
       NLL+=log(evaluated);
    }
    return  -NLL;
};

//Method to get the normalization value
double myGauss::normValue(){
     double  output = paramValue[0]*sqrtTwoPi;
     for(int i=1; i<dimensions; i++){
         output *= paramValue[i]*sqrtTwoPi;
     };
     return output;
};


//Method to generate data distributed like the gaussian
void myGauss::generateData(int length, double *  p){
    srand(time(NULL));
    for(int i=0; i<length*dimensions; i++){
        double rand1 = ((double)rand())/RAND_MAX;
        double rand2 = ((double)rand())/RAND_MAX;
        //Now doing the box muller method
        p[i] = sqrt(-2.0*log(rand1))*cos(2*PI*rand2)*paramValue[i%dimensions];
    };
};

//Method to get number of dimensions for  function
int myGauss::getDimensions(){
    return dimensions;
};
