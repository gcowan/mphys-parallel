//Class to do negative log likelyhood
#ifndef NEGATIVELOG_H
#define NEGATIVELOG_H
class negativeLog{
    public:
        //Method get the negative log likelihood from the dataSet  
        double evaluateDataSet(double * data, int numData, double sigmaEst);
};
#endif

