#include <cilk/cilk.h>
#ifndef MYFUNC_H
#define MYFUNC_H
class _Cilk_shared myFunc{
    protected: 
        double * paramValue;
        int dimensions;
    public:
         virtual void  generateData(int outputLength,double * p)=0;
         virtual double  evaluate( double *  dataSet, int dataLength, int threads)=0;
         virtual void setParameter(int paraIndex, double para)=0;
         virtual double getParameter(int paraIndex)=0;
         virtual double normValue()=0;
         virtual int getDimensions()=0;
};
#endif
