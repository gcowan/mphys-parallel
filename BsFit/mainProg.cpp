#define PI 3.141592653589793238462
#define  FREE  alloc_if(0) free_if(1)
#define  REUSE  alloc_if(0) free_if(0)
#define  ALLOC  alloc_if(1) free_if(0)

#pragma offload_attribute(push,target(mic))
double * params;
double * dataSet;
int dataLength=235764;
#include <omp.h>
#include <stdio.h>
#pragma offload_attribute(pop)
//Number of parameters
int parDim = 2;
#include "fit.h"
//Minuit interface
#include <mkl_vsl.h>
#include "TMinuit.h"
void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag);



int main(){
    params = (double *)malloc(sizeof(double)*parDim);
    for(int i=0; i<parDim; i++){
        params[i] = 0;
    }
    FILE * input = fopen("BsToJpsiPhi_LHCb_data.txt","r");
    dataSet = (double *) malloc(sizeof(double)*dataLength*dimensions);
    for(int i=0; i<dataLength; i++){
       double tmp = 0;
       //Parsing masses ignoring lifetimes
       fscanf(input,"%lf",&tmp);
       dataSet[i]=tmp;
       //Now scanning lifetime and doing nothing with it
       fscanf(input,"%lf",&tmp);
    }
    printf("Finished scanning  data element 1 is %f should be 5485.2739  \n",dataSet[0]);
    printf("data element 2 is %f should be 5422.9331\n",dataSet[1]);
    //Allocating dataSet on coprocessor
        #pragma offload_transfer target(mic:0) in(dataSet : length(dataLength*dimensions)  ALLOC) 
        #pragma offload_transfer target(mic:0) in( dataLength : ALLOC)
    double time = omp_get_wtime();
    TMinuit min(parDim);
    min.SetFCN(funcn);
    min.SetErrorDef(0.5);
    min.SetPrintLevel(1); 
    min.DefineParameter(0,"Mass", 5200, 100,5000,5500);
    min.DefineParameter(1," Width of Signal",50,5,0,200);
     // min.DefineParameter(2,"Background frac",5,1,0,10);
     // min.DefineParameter(3,"Background A",0,1,-5,5);
/*  min.DefineParameter(5,"Sigma 5",0.7,0.1,0.5,1);
    min.DefineParameter(6,"Sigma 6",0.7,0.1,0.5,1);
    min.DefineParameter(7,"Sigma 7",0.7,0.1,0.5,1);
    min.DefineParameter(8,"Sigma 8",0.7,0.1,0.5,1);
    min.DefineParameter(9,"Sigma 9",0.7,0.1,0.5,1);
    */
    printf("Minuit Start\n");
    min.Migrad();
    
    time = omp_get_wtime()-time;
    printf("That fit took %f seconds \n",time);
    fclose(input);
    //Freeing dataSet
        #pragma offload_transfer target(mic:0) nocopy(dataSet : FREE)
    return 0;
}




void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    params[0] = par[0];
    params[1] = par[1];
     // params[2] = par[2];
     // params[3] = par[3];
    double result = 0;
     // printf("Func call\n");
     // double time = omp_get_wtime();
          #pragma offload target(mic:0) in(params : length(parDim)) out(result) nocopy(dataSet: length(0) REUSE)
    {
            result = evaluateDataSet(dataSet, dataLength,240,params);
    } 
     // time = omp_get_wtime()-time;
     // printf("Took %f seconds\n",time);
    f = result;
       // printf("result=%f\n",result);
       // printf("par1=%f \n",par[0]);
};



