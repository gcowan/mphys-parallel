#include "myFunc.h"
#include "myGauss.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <cilk/cilk.h>
#pragma offload_attribute(push,_Cilk_shared)
#include <offload.h>
#pragma offload_attribute(pop)
#include "TMinuit.h"


using namespace std;
_Cilk_shared  myFunc * _Cilk_shared  funcp;
_Cilk_shared const int length = 1000000;
_Cilk_shared  double  *  _Cilk_shared  data;  
_Cilk_shared int threads=12;
int dimensions = 10;
_Cilk_shared myGauss gauss(dimensions);
_Cilk_shared double * _Cilk_shared limits;
 FILE * output;

//try to trick the compiler
double _Cilk_shared wrapperFunction(){
    return funcp->evaluate(data,length,threads);
};

void _Cilk_shared funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    int dim = funcp->getDimensions();
    for(int i=0; i<dim; i++){
        funcp->setParameter(i,par[i]);
    };
    double result =     wrapperFunction();
    f = result;
       // printf("result=%f\n",result);
       // printf("par1=%f \n",par[0]);
};

double _Cilk_shared wrapperMC(_Cilk_shared double * lims){
    return funcp->integrateVegas(lims,240);
}


int main(){
    funcp = &gauss;
    for(int i=0; i<dimensions; i++){
        funcp->setParameter(i,0.783);
    };
      data =  (_Cilk_shared double * ) _Offload_shared_aligned_malloc(sizeof(double)*length*dimensions,64);
      funcp->generateData(length,data);
    printf("Finished generating data\n");
    double time = omp_get_wtime();
    TMinuit min(dimensions);
    min.SetFCN(funcn);
    min.SetErrorDef(0.5);
    min.SetPrintLevel(1); 
    min.DefineParameter(0,"Sigma", 0.7, 0.1,0.5,1.0);
    min.DefineParameter(1,"Sigma 1",0.7,0.1,0.5,1.0);
    min.DefineParameter(2,"Sigma 2",0.7,0.1,0.5,1);
    min.DefineParameter(3,"Sigma 3",0.7,0.1,0.5,1);
    min.DefineParameter(4,"Sigma 4",0.7,0.1,0.5,1);
    min.DefineParameter(5,"Sigma 5",0.7,0.1,0.5,1);
    min.DefineParameter(6,"Sigma 6",0.7,0.1,0.5,1);
    min.DefineParameter(7,"Sigma 7",0.7,0.1,0.5,1);
    min.DefineParameter(8,"Sigma 8",0.7,0.1,0.5,1);
    min.DefineParameter(9,"Sigma 9",0.7,0.1,0.5,1);
    min.Migrad();
    time = omp_get_wtime()-time;
    printf("That fit took %f seconds \n",time);
    return 0;
};
