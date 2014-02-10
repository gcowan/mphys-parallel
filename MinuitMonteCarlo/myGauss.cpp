#include "myGauss.h"
#include <math.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <mkl_vsl.h>
#define PI 3.141592653589793238462
//Number of points to sample in monte carlo integration
#define nsamples 100000

//Generate Gaussian with set of parameters
myGauss::myGauss(int numSigmas,  double *  paramValues){
    paramValue = paramValues;
    dimensions = numSigmas;
    //CHANGE SEED WHEN YOU KNOW IT WORKS
    vslNewStream(&stream, VSL_BRNG_SFMT19937,3864368);
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
    //CHANGE SEED WHEN YOU KNOW IT WORKS
    vslNewStream(&stream, VSL_BRNG_SFMT19937,3864368);
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
     srand(time(NULL));
     double integral=0;
     double eval;
     int lower = -1;
     int higher = 1;
     //Doing integration in parallel
     #pragma omp parallel for default(none) private(eval) shared(lower,higher) reduction(+:integral)
     for(int i=0; i<nsamples; i++){
         //Evaluation of function with random values of X
         eval = 1;
         for(int j=0; j<dimensions; j++){
             double randNum = randDouble(lower,higher);
             eval*=exp(-randNum*randNum*oneOverTwoSigsSq[j]);
         }

         integral+=eval;
     };
     integral /= nsamples;
     integral*=(higher-lower);

     double  output = paramValue[0]*sqrtTwoPi;
     for(int i=1; i<dimensions; i++){
         output *= paramValue[i]*sqrtTwoPi;
     };
     printf("Integral is:%f\n",integral);
     printf("Maths prediction=%f\n",output);
     return integral;
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

//Method to return random doubles generated between limits
double myGauss::randDouble(int lower, int higher){ 
    double val;
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,&val,lower,higher);
    return val;
};

//Integration using the vegas algorithm
double myGauss::integrateVegas(double * limits){
    //How many iterations to perform
    int iterations = 10;
    //How many points to sample in total
    int samples = 10000000;
    //How many intervals for each dimension
    int intervals = 10;
    //How many subIntervals
    int subIntervals = 1000;
    //Parameter alpha controls convergence rate
    double alpha = 0.2;
    //Number of boxes
    int numBoxes = intervals;
    for(int i=1; i<dimensions; i++){
        numBoxes *= intervals;
    }
    //Arrays to store integral and uncertainty for each iteration
    double * integral = (double *)_Offload_shared_aligned_malloc(sizeof(double)*iterations,64);
    double * sigmas = (double *)_Offload_shared_aligned_malloc(sizeof(double)*iterations,64);
    //Points per each box
    int pointsPerBox = samples/numBoxes;
    //Array storing the box limits (stores x limits then y limits and so on) intervals+1 to store all limits
    double * boxLimits = (double *)_Offload_shared_aligned_malloc(sizeof(double)*(intervals+1)*dimensions,64);
    //Array to store average function values for each box
    double * heights = (double *)_Offload_shared_aligned_malloc(sizeof(double)*dimensions*intervals,64);
    //Array storing values of m
    double * mValues = (double *)_Offload_shared_aligned_malloc(sizeof(double)*intervals,64);
    //Array storing widths of sub boxes
    double * subWidths = (double *) _Offload_shared_aligned_malloc(sizeof(double)*intervals,64);
    //Getting initial limits for the boxes 
    for(int i=0; i<dimensions; i++){
        double boxWidth = (limits[i+1]-limits[i])/intervals;
        //0th iteration
        boxLimits[i*(intervals+1)] = limits[i];
        for(int j=1; j<=intervals; j++){
            int x = (i*(intervals+1))+j;
            boxLimits[x] =  boxLimits[x-1]+boxWidth;
        }
    };
    //Pointer to store random generated  numbers
    double * randomNums = (double *)_Offload_shared_aligned_malloc(sizeof(double)*dimensions,64);
    int * binNums = (int *)_Offload_shared_aligned_malloc(sizeof(int)*dimensions,64);
    for(int iter=0; iter<iterations; iter++){ 
        //Performing  iterations
        for(int i=0; i<samples; i++){
            viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,dimensions,binNums,0,intervals);
            //Getting samples from bins
            for(int j=0; j<dimensions; j++){
                int x = ((intervals+1)*j)+binNums[j];
                randomNums[j] = randDouble(boxLimits[x],boxLimits[x+1]);
            }
            //Performing evaluation of function and adding it to the total integral
            double eval = evaluate(randomNums);
            integral[iter] += eval;
            sigmas[iter] += (eval*eval);
            //Calculating the values of f for bin resising
            for(int j=0; j<dimensions; j++){
                int x = binNums[j]+(j*intervals);
                //May need to initialize heights
                heights[x] += eval;
            }

        } 
        //Calculating the values of sigma and the integral
        integral[iter] /= samples;
        sigmas[iter] /= samples;
        sigmas[iter] -= (integral[iter]*integral[iter]);
        sigmas[iter] /= (samples-1);

        //Readjusting the box widths based on the heights
        //Creating array to store values of m and their sum 
        int totalM=0; 
        //Doing for each dimension seperately
        //May be a bug in this for loop TEST
        for(int i=0; i<dimensions; i++){
            double sum = 0;
            //Getting the sum of f*delta x
            for(int j=0; j<intervals; j++){
                int x = (i*(intervals))+j ;
                //May be bug with these indicies
                sum += heights[x]*(boxLimits[x+1+i]-boxLimits[x+i]);
            }
            //Performing the rescaling 
            for(int j=0; j<intervals; j++){
                int x = (i*(intervals))+j;
                double value = heights[x]*(boxLimits[x+1+i]-boxLimits[x+i]);
                mValues[j] = ceil(subIntervals*pow((value-1)*(1.0/log(value)),alpha));
                subWidths[j] = boxLimits[x+1+i]-boxLimits[x+i]/mValues[j];
                totalM += mValues[j];
            }
            int mPerInterval = totalM/intervals;
            int mValueIterator = 0;
            //Adjusting the intervals going from 1 to less than intervals to keep the edges at the limits
            for(int j=1; j<intervals; j++){
                double width = 0;
                for(int y=0; y<mPerInterval; y++){
                    width += subWidths[mValueIterator];
                    mValues[mValueIterator]--;
                    if(mValues[mValueIterator]==0){
                        mValueIterator++;
                    }
                }
                //NEED TO SET BOX LIMITS NOW  
                int x = j+(i*(intervals+1));
                boxLimits[x] = boxLimits[x-1]+width;    
            }
            //Setting mvalues etc. (reseting memory allocated before the dimensions loop to 0)
            totalM = 0;
            for(int k=0; k<intervals; k++){
                subWidths[k] = 0;
                mValues[k] = 0;

            }
        }
        //Setting heights to zero for next iteration
        for(int i=0; i<intervals*dimensions; i++ ){
            heights[i] = 0;
        }
    }

    //All iterations done 
    //Free stuff
    _Offload_shared_aligned_free(subWidths);
    _Offload_shared_aligned_free(mValues);
    _Offload_shared_aligned_free(binNums);
    _Offload_shared_aligned_free(randomNums);
    _Offload_shared_aligned_free(boxLimits);
    //Calculating the final value of the integral
    double denom = 0;
    double numerator =0;
    for(int i=0; i<iterations; i++){
        numerator += integral[i]*(integral[i]*integral[i]/(sigmas[i]*sigmas[i]));
        denom += (integral[i]*integral[i]/(sigmas[i]*sigmas[i]));
    }
    double output  = numerator/denom;
    //Calculating value of x^2 to check if result can be trusted
    double chisq = 0;
    for(int i=0; i<iterations; i++){
       chisq += ((integral[i]-output)*(integral[i]-output)/(output*output))*(integral[i]*integral[i]/(sigmas[i]*sigmas[i]));
    }
    if(chisq>iterations){
        printf("Chisq value is %f, it should be not much greater than %d (iterations-1)\n",chisq,iterations-1);
    }
    _Offload_shared_aligned_free(integral);
    _Offload_shared_aligned_free(sigmas);
    return output;
    
}

