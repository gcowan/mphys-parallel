#include "negativeLog.h"
#include "gaussianFunc.h"
#include <cmath>
#include <omp.h>

double negativeLog::evaluateDataSet(double * data, int dataLength, double sigmaVal){
    gaussianFunc func = gaussianFunc(sigmaVal);
    double normalization = func.integral();
    //summing up all the logs
    int numThreads = 8;
    double total [numThreads];
    #pragma omp parallel  default(none) shared(total,numThreads,normalization,dataLength,data,func)
    {
        int threadNum = omp_get_thread_num();
        int lowerBound = (dataLength/numThreads)*(threadNum);
        int upperBound =  (dataLength/numThreads)*(threadNum+1);

        for(int i =lowerBound; i<upperBound; i++){
            double val = func.evaluate(data[i]);
            total[threadNum] += log(val/normalization);
        }
        //Making sure we iterate upto the last data element
        if(threadNum==numThreads-1 && upperBound<dataLength){
            for(int i=upperBound; i<dataLength; i++){
                double val = func.evaluate(data[i]);
                total[threadNum] += log(val/normalization);
            }
        }

    }
    //Now summing up the totals produced by each thread
    double trueTotal = 0;
    for(int i=0; i<numThreads; i++){
        trueTotal += total[i];
    }
    return -trueTotal;
};
