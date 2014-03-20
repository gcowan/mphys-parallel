#define PI 3.141592653589793238462
#define dimensions 1
#define  FREE  alloc_if(0) free_if(1)
#define  REUSE  alloc_if(0) free_if(0)
#define  ALLOC  alloc_if(1) free_if(0)

//Headers methods run on both sides
#pragma offload_attribute(push,target(mic))
double evaluate(double * value, double * oneOverTwoSigsSq);
double evaluateDataSet(double * dataSet, int dataLength,int threads,double * params);
inline double integrateVegas(double * limits, int threads, double * params);
#pragma offload_attribute(pop)
