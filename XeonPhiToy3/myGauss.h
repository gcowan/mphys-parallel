#include "myFunc.h"
#include <cilk/cilk.h>
#ifndef MYGAUSS_H
#define MYGAUSS_H
class _Cilk_shared myGauss: public myFunc {
    private:
        double oneOverTwoSigSq;
        double sqrtTwoPi;
    public:
        myGauss(double sigma);
        myGauss();
        ~myGauss();
        double getParameter();
        void  setParameter(double newSigma);
        void generateData(int outputLength , double * p);
        double evaluate(double xValue);
        double evaluate(_Cilk_shared double * dataSet, int dataLength, int threads);
        double normValue();

};
#endif
