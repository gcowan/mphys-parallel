#include "myFunc.h"
#include <cilk/cilk.h>
#ifndef MYEXP_H
#define MYEXP_H
class _Cilk_shared myExp : public  myFunc{
    public:
        myExp(double expValue);
        myExp();
        ~myExp();
        double getParameter();
        void setParameter(double newExp);
        void generateData(int dataLength, double * p);
        double evaluate(_Cilk_shared double *  data, int dataLength, int threads);
        double normValue();


};
#endif
