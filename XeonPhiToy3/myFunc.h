#include <cilk/cilk.h>
#ifndef MYFUNC_H
#define MYFUNC_H
#pragma offload_attribute (push,_Cilk_shared)
class myFunc{
    protected: 
        double paramValue;
    public:
        myFunc() {paramValue = 1; };
        myFunc(double param){
            paramValue = param;
        };
        ~myFunc(){};
        virtual double * generateData(int outputLength)=0;
        virtual double  evaluate(double * dataSet, int dataLength)=0;
        virtual void setParameter(double para){paramValue = para;};
        virtual double getParameter(){return paramValue;};
        virtual double normValue()=0;
};
#pragma offload_attribute(pop)
#endif
