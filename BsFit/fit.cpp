#define PI 3.141592653589793238462
#define  FREE  alloc_if(0) free_if(1)
#define  REUSE  alloc_if(0) free_if(0)
#define  ALLOC  alloc_if(1) free_if(0)


//Headers methods run on both sides
#pragma offload_attribute(push,target(mic))
#include "fit.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mkl_vsl.h>
#pragma offload_attribute(pop)


//Method to evaluate function for one set of x values
inline double evaluate(double * value, double * params){
    //Gaussian Signal
    double expon = value[0]-params[0];

    double output = exp(-0.5*expon*expon/(params[1]*params[1]));
    //Background Noise polynomial
     // output += params[2]*(1+params[3]*value[0]);
    // printf("output=%f value[0]=%f params[0]=%f sigma=%f\n",output,value[0],params[0],params[1]);
    return output;
};


//Method to evaluate entire dataset using a set number of threads
double evaluateDataSet( double * dataSet, int dataLength, int threads, double * params){
    double NLL = 0;
    double * limits = (double *) malloc(sizeof(double)*dimensions*2);
    for(int i=0; i<dimensions; i++){
        limits[2*i] = 0.0;
        limits[2*i+1] = 10000.0;
    }
    double normTime = omp_get_wtime();
    double norm = 1.0/integrateVegas(limits,threads,params);
    normTime = omp_get_wtime()-normTime;
    printf("Time to normalize = %f\n",normTime);
    // double norm = 1.0/normValue(params);
    omp_set_num_threads(threads);
    /*int off = _Offload_get_device_number();
      if(off==0)
      printf("On mic\n");
      if(off==-1)
      printf("On host\n");
      */
    double evalTime = omp_get_wtime();
#pragma omp parallel for default(none) shared(dataSet,dataLength,norm, params) reduction(+:NLL)
    for(int i=0; i<dataLength*dimensions; i+=dimensions){
        double evaluated = evaluate(dataSet+i,params);
        NLL+=log(evaluated);
    }
    evalTime = omp_get_wtime()-evalTime;
    printf("Time to evaluate %f\n",evalTime);
    free(limits);
    return  -NLL;
};








inline double integrateVegas(double * limits , int threads, double * params){
    //Setting the number of threads
    omp_set_num_threads(threads);
    //How many iterations to perform
    int iterations =15;
    //Which iteration to start sampling more
    int switchIteration = 7;
    //How many points to sample in total
    int samples = 100000;
    //How many points to sample after grid set up
    int samplesAfter = 5000000;
    //How many intervals for each dimension
    int intervals = 10;
    //How many subIntervals
    int subIntervals = 1000;
    //Parameter alpha controls convergence rate
    double alpha = 0.5;
    int seed = 92834721;
    //double to store volume integrated over
    double volume = 1.0;
    for(int i=0; i<dimensions; i++){
        volume*= (limits[(2*i)+1]-limits[2*i]);
    };
    //Number of boxes
    int numBoxes = intervals;
    for(int i=1; i<dimensions; i++){
        numBoxes *= intervals;
    }
    //CHANGE SEED WHEN YOU KNOW IT WORKS
    //Setting up one random number stream for each thread
    VSLStreamStatePtr * streams; 
    streams = ( VSLStreamStatePtr * )malloc(sizeof(VSLStreamStatePtr)*threads);
    for(int i=0; i<threads; i++){
        vslNewStream(&streams[i], VSL_BRNG_MT2203+i,seed);
    }
    //Arrays to store integral and uncertainty for each iteration
    double * integral = (double *)malloc(sizeof(double)*iterations);
    double * sigmas = (double *)malloc(sizeof(double)*iterations);
    for(int i=0; i<iterations; i++){
        integral[i] = 0;
        sigmas[i] = 0;
    }
    //Points per each box
    int pointsPerBox = samples/numBoxes;
    //Array storing the box limits (stores x limits then y limits and so on) intervals+1 to store all limits
    double * boxLimits = (double *)malloc(sizeof(double)*(intervals+1)*dimensions);
    //Array to store average function values for each box
    double * heights = (double *)malloc(sizeof(double)*dimensions*intervals);
    //Array storing values of m
    double * mValues = (double *)malloc(sizeof(double)*intervals);
    //Array storing widths of sub boxes
    double * subWidths = (double *) malloc(sizeof(double)*intervals);
    //Getting initial limits for the boxes 
    for(int i=0; i<dimensions; i++){
        double boxWidth = (limits[(2*i)+1]-limits[2*i])/intervals;
        //0th iteration
        boxLimits[i*(intervals+1)] = limits[2*i];
        for(int j=1; j<=intervals; j++){
            int x = (i*(intervals+1))+j;
            boxLimits[x] =  boxLimits[x-1]+boxWidth;
        }
    };
    //Pointer to store random generated  numbers
    __declspec(align(64)) double  randomNums [dimensions];
    __declspec(align(64))  int  binNums [dimensions];
    //Double to store p(x) denominator for monte carlo
    double prob;
    //Values to store integral and sigma for each thread so they can be reduced in OpenMp
    double integralTemp;
    double sigmaTemp;
    __declspec(align(64)) double heightsTemp[dimensions*intervals];
    int threadNum;
    double timeMC; 
#pragma omp parallel  default(none) private(sigmaTemp,integralTemp,binNums,randomNums,prob,threadNum,heightsTemp) shared(iterations,subIntervals,alpha,mValues,subWidths,streams,samples,boxLimits,intervals, integral, sigmas, heights, threads, volume, samplesAfter, switchIteration, params,timeMC) 
    {

        for(int iter=0; iter<iterations; iter++){ 

#pragma omp single
            {
                timeMC= omp_get_wtime();
            }
            //Stepping up to more samples when grid calibrated
            if(iter==switchIteration){
                samples = samplesAfter;
            }
            //Performing  iterations
            for(int i=0; i<dimensions*intervals; i++){
                heightsTemp[i] = 0;
            }

            integralTemp = 0; 
            sigmaTemp = 0;
            //Getting chunk sizes for each thread
            threadNum = omp_get_thread_num();
            int seg = ceil((double)samples/threads);
            int lower = seg*threadNum;
            int upper = seg*(threadNum+1);
            if(upper > samples){
                upper = samples;
            };
            //Spliting monte carlo up
            for(int i=0; i<seg; i++){
                prob = 1;
                //Randomly choosing bins to sample from
                viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,streams[threadNum],dimensions,binNums,0,intervals);
                vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,streams[threadNum],dimensions,randomNums,0,1);
                //Getting samples from bins
                for(int j=0; j<dimensions; j++){
                    int x = ((intervals+1)*j)+binNums[j];
                    randomNums[j] *= (boxLimits[x+1]-boxLimits[x]);
                    randomNums[j] += boxLimits[x];
                    prob *= 1.0/(intervals*(boxLimits[x+1]-boxLimits[x]));
                }
                //Performing evaluation of function and adding it to the total integral
                double eval = evaluate(randomNums,params);
                integralTemp += eval/prob;
                sigmaTemp += (eval*eval)/(prob*prob);
                //Calculating the values of f for bin resising
                for(int j=0; j<dimensions; j++){
                    int x = binNums[j]+(j*intervals);
                    //May need to initialize heights
                    // #pragma omp atomic
                    // printf("heightsTemp before=%f\n",heightsTemp[x]);
                    heightsTemp[x] += eval;
                    // printf("heightsTemp=%f x=%d eval=%f thread=%d\n",heightsTemp[x],x,eval,omp_get_thread_num());
                }

            } 
#pragma omp critical
            {
                integral[iter] += integralTemp;
                sigmas[iter] += sigmaTemp;
                for(int k=0; k<dimensions*intervals; k++){
                    heights[k] += heightsTemp[k];
                }
            }
#pragma omp barrier
            // double resizeTime = omp_get_wtime();
#pragma omp single
            {
                //Calculating the values of sigma and the integral
                integral[iter] /= samples;
                sigmas[iter] /= samples;
                sigmas[iter] -= (integral[iter]*integral[iter]);
                sigmas[iter] /= (samples-1);
                 // printf("integral=%f\n",integral[iter]);

                //Readjusting the box widths based on the heights
                //Creating array to store values of m and their sum 
                int totalM=0; 
                //Doing for each dimension seperately
                double loop2Time = 0;
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
                        double value = heights[x]*(boxLimits[x+1+i]-boxLimits[x+i])/sum;
                        mValues[j] = ceil(subIntervals*pow((value-1)*(1.0/log10(value)),alpha));
                        subWidths[j] = (boxLimits[x+1+i]-boxLimits[x+i])/mValues[j];
                        totalM += mValues[j];
                    }
                    int mPerInterval = totalM/intervals;
                    int mValueIterator = 0;
                    //Adjusting the intervals going from 1 to less than intervals to keep the edges at the limits
                    double loo2Time = omp_get_wtime();
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
                    loo2Time = omp_get_wtime()-loo2Time;
                    //Setting mvalues etc. (reseting memory allocated before the dimensions loop to 0)
                    totalM = 0;
                    for(int k=0; k<intervals; k++){
                        subWidths[k] = 0;
                        mValues[k] = 0;

                    }
                    loop2Time += loo2Time;
                }

                //Setting heights to zero for next iteration
                for(int i=0; i<intervals*dimensions; i++ ){
                    heights[i] = 0;
                }

                 // printf("Single time=%f  loop2Time=%f\n",singleTime,loop2Time);
            }
            // resizeTime = omp_get_wtime()-resizeTime;
            // printf("Resizing the boxes too %f seconds\n",resizeTime);


        }
    }
    //All iterations done 
    //Free stuff
    
    free(subWidths);
    free(mValues);
    free(boxLimits);
    free(streams);
    free(heights);
   
    //Calculating the final value of the integral
    double denom = 0;
    double numerator =0;
    for(int i=7; i<iterations; i++){
        numerator += integral[i]*((integral[i]*integral[i])/(sigmas[i]*sigmas[i]));
        denom += ((integral[i]*integral[i])/(sigmas[i]*sigmas[i]));
         // printf("integral=%f sigma=%f\n",integral[i],sigmas[i]);
    }
    double output  = numerator/denom;
    //Calculating value of x^2 to check if result can be trusted
    double chisq = 0;
    for(int i=0; i<iterations; i++){
       chisq += (((integral[i]-output)*(integral[i]-output))/(sigmas[i]*sigmas[i]));
    }
     // if(chisq>iterations){
         // printf("Chisq value is %f, it should be not much greater than %d (iterations-1) Integral:%f \n",chisq,iterations-1,output);
     // }
      free(integral);
      free(sigmas);
    return output;
    
}

