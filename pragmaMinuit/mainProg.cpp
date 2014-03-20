#define PI 3.141592653589793238462
#define  FREE  alloc_if(0) free_if(1)
#define  REUSE  alloc_if(0) free_if(0)
#define  ALLOC  alloc_if(1) free_if(0)

#pragma offload_attribute(push,target(mic))
double * params;
double * dataSet;
int dataLength=1000000;
int threads;
double timeIn=0;
#include <omp.h>
#include <stdio.h>
#pragma offload_attribute(pop)

#include "fit.h"
//Minuit interface
#include <mkl_vsl.h>
#include "TMinuit.h"
void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag);


void generateData(int length, double *  p, double * paramValue){
    VSLStreamStatePtr stream;
    vslNewStream( &stream, VSL_BRNG_MT19937, 239107 );
    for(int i=0; i<length*dimensions; i++){
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,1,(p+i),0,paramValue[i%dimensions]);
    };
};


int main(){
    params = (double *)malloc(sizeof(double)*dimensions);
    for(int i=0; i<dimensions; i++){
        params[i] = 0.783;
    }
    dataSet = (double *) malloc(sizeof(double)*dataLength*dimensions);
    generateData(dataLength,dataSet,params);
    printf("Finished generating data\n");
    FILE * output = fopen("hostMC.txt","w");
    double mallocTime = 0;
    double singleTime = 0;
    double phiTime = 0;
    //Allocating dataSet on coprocessor
      // #pragma offload_transfer target(mic:0) in(dataSet : length(dataLength*dimensions)  ALLOC) 
       // #pragma offload_transfer target(mic:0) in( dataLength : ALLOC)
    for(threads=1; threads<=12; threads+=1){
/*        TMinuit min(dimensions);
        min.SetFCN(funcn);
        min.SetErrorDef(0.5);
        min.SetPrintLevel(1); 
        min.DefineParameter(0,"Sigma", 0.7, 0.1,0.5,1.0);
        min.DefineParameter(1,"Sigma 1",0.7,0.1,0.5,1.0);
        min.DefineParameter(2,"Sigma 2",0.7,0.1,0.5,1);
        min.DefineParameter(3,"Sigma 3",0.7,0.1,0.5,1);
        min.DefineParameter(4,"Sigma 4",0.7,0.1,0.5,1);
        double overallTime = omp_get_wtime();
        min.Migrad();
        overallTime = omp_get_wtime()-overallTime;
        int numCalls =  min.fNfcn;
        fprintf(output,"%d %f %f %f %f\n",threads,overallTime,timeIn,overallTime/numCalls,timeIn/numCalls);
        timeIn = 0;
        */
        double * limits = (double *) malloc(sizeof(double)*dimensions*2);
        for(int i=0; i<dimensions; i++){
            limits[2*i] = -5.0;
            limits[2*i+1] = 5.0;
        }
        double result;
         // #pragma offload_transfer target(mic:0) in( limits : length(2*dimensions) ALLOC)
 // #pragma offload_transfer target(mic:0) in(params : length(dimensions) ALLOC)
        double totalTime = omp_get_wtime();
       
 // #pragma offload target(mic:0) out(result) inout(mallocTime,singleTime,phiTime) nocopy(limits,params : REUSE)
        {
            result = integrateVegas(limits,threads,params,mallocTime,singleTime,phiTime);
        }
        totalTime = omp_get_wtime()-totalTime;
        fprintf(output,"%d %f %f %f %f\n",threads,totalTime,phiTime,mallocTime,singleTime);
        printf("result=%f normval=%f\n",result,normValue(params));
        mallocTime =0;
        singleTime =0;
        
    } 
    fclose(output);
    //Freeing dataSet
       // #pragma offload_transfer target(mic:0) nocopy(dataSet : FREE)
    return 0;
}




void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    for(int i=0; i<dimensions; i++){
        params[i] = 1.0/(2*par[i]*par[i]);
    }; 
    double result = 0;
         // #pragma offload target(mic:0) in(params : length(dimensions)) out(result) nocopy(dataSet: length(0) REUSE) inout(timeIn)
    {
            result = evaluateDataSet(dataSet, dataLength,threads,params,timeIn);
    } 
    f = result;
       // printf("result=%f\n",result);
       // printf("par1=%f \n",par[0]);
};



