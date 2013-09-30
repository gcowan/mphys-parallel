#include "negativeLog.h"
#include "guassianFunc.h"
#include <cmath>

double negativeLog::evaluateDataSet(double * data, int dataLength, double sigmaVal){
    gaussianFunc func = gaussianFunc(sigmaVal);
    double normalization = func.integral();
    //summing up all the logs
    double total = 0;
    for(int i =0; i<dataLength; i++){
        double val = func.evaluate(data[i]);
        total+= log(val/normalization);
    }
    return -total;
};
