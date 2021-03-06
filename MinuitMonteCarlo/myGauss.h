#include "myFunc.h"
#include <cilk/cilk.h>
#ifndef MYGAUSS_H
#define MYGAUSS_H
class _Cilk_shared myGauss: public myFunc {
    private:
        double * oneOverTwoSigsSq;
        double sqrtTwoPi;

    public:
        double avgTime;
        myGauss(int numSigmas , double *   sigmas);
        myGauss(int numSigmas);
        ~myGauss();
        double getParameter(int index);
        void  setParameter(int index,double newSigma);
        void generateData(int outputLength , double * p);
        double evaluate(double * values);
        double evaluate( double * dataSet, int dataLength, int threads);
        double normValue();
        int getDimensions();
        void generateRandNums(double *, int *, int,int);
        double integrateVegas( double * , int threads );


};
#endif
