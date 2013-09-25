#include "numArray.h"
#include "numArrayWrapper.h"
#include <omp.h>
#include <cstdio>

double numArrayWrapper::sqrtAverageFast(numArray data){
    double average = 0;
    data.sqrtParallelArray();
    int numThreads = omp_get_num_threads();
    int threadNum;
    int length = data.getLength();
    int i;
    printf("numThreads=%d\n",numThreads);
    #pragma omp parallel shared(average,numThreads,length,data) private(threadNum,i) default(none)
    { 
    	//getting the thread id
	threadNum = omp_get_thread_num();
	int startValue = (length/numThreads)*threadNum;
	int endValue = (length/numThreads)*(threadNum+1)-1;
	printf("startValue=%d   endValue=%d  threadNum=%d\n", startValue,endValue, threadNum);
        for( i=startValue; i<endValue; i++){
            average+=data.getElement(i);
	    
	    if(i==length){
	         printf("Got to the end of the data on thread %d\n",threadNum);
	    }
        }
    }
    return average/data.getLength();
    
};

double numArrayWrapper::sqrtAverageSlow(numArray data){
    double average = 0;
    data.sqrtArray();
    for(int i=0; i<data.getLength(); i++){
        average+=data.getElement(i);
    }
    return average/data.getLength(); 
    
};
