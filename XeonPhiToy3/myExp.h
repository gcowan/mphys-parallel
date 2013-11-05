#include "myFunc.h"
#include <cilk/cilk.h>
#ifndef MYEXP_H
#define MYEXP_H
#pragma offload_attribute (push,_Cilk_shared) 
class myExp : public  myFunc{
    public:
        myExp(double expValue);
        myExp();
        ~myExp();
        double getParameter();
        void setParameter(double newExp);
        double * generateData(int dataLength);
        double evaluate(double * data, int dataLength);
        double normValue();


};
#pragma offload_attribute (pop)
#endif
