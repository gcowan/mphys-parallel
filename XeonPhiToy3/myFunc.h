#include <cilk/cilk.h>
#ifndef MYFUNC_H
#define MYFUNC_H
class _Cilk_shared myFunc{
    protected: 
        double paramValue;
    public:
         myFunc() {paramValue = 1; };
         myFunc(double param){
            paramValue = param;
        };
        ~myFunc(){};
         virtual void  generateData(int outputLength,double * p)=0;
         virtual double  evaluate(_Cilk_shared double *  dataSet, int dataLength, int threads)=0;
         virtual void setParameter(double para){paramValue = para;};
         virtual double getParameter(){return paramValue;};
         virtual double normValue()=0;
};
#endif
