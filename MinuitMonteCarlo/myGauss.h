#include "myFunc.h"
#include <cilk/cilk.h>
#include <mkl_vsl.h>
#ifndef MYGAUSS_H
#define MYGAUSS_H
class _Cilk_shared myGauss: public myFunc {
    private:
        double * oneOverTwoSigsSq;
        double sqrtTwoPi;
        VSLStreamStatePtr * streams;

    public:
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
        double randDouble(double,double, int thread);
        double integrateVegas( double * , int threads );


};
#endif
