#include "myFunc.h"
#include <cilk/cilk.h>
#ifndef MYGAUSS_H
#define MYGAUSS_H

#pragma offload_attribute (push,_Cilk_shared)
class myGauss: public myFunc {
    private:
        double oneOverTwoSigSq;
        double sqrtTwoPi;
    public:
        myGauss(double sigma);
        myGauss();
        ~myGauss();
        double getParameter();
        void  setParameter(double newSigma);
        double * generateData(int outputLength);
        double evaluate(double xValue);
        double evaluate(double * dataSet, int dataLength);
        double normValue();

};

#pragma offload_attribute(pop)
#endif
