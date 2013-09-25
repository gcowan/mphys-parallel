#include "numArray.h"
#include "numArrayWrapper.h"
#include <omp.h>
#include <cstdio>

double numArrayWrapper::sqrtAverageLogFast(numArray data){
    data.sqrtParallelArray();
    data.logParallelArray();
    int numThreads = 8;
    double averageArray[numThreads];
    int threadNum;
    int length = data.getLength();
    int i;
    
    for(int i=0; i<numThreads; i++){
       averageArray[i] = 0; 
    }


    //printf("numThreads=%d\n",numThreads);
    #pragma omp parallel shared(numThreads,length,data,averageArray) private(threadNum,i) default(none)
    { 
    	//getting the thread id
	threadNum = omp_get_thread_num();
	//Splitting up the array between the threads
	int startValue = (length/numThreads)*threadNum;
	int endValue = (length/numThreads)*(threadNum+1);
	int range = endValue-startValue;
	//printf("threadnum %d has range = %d",threadNum,range);
        for( i=startValue; i<endValue; i++){
            averageArray[threadNum]+=data.getElement(i);
        }
	//Doing a check for the last variables to make sure they are included
	if(threadNum==numThreads-1 && endValue<length){
	    //getting difference between end readup to by last thread and end of data
	    for(int i=endValue; i<length; i++){
	        averageArray[threadNum]+=data.getElement(i);
	    }
	}
    }
    
    double average = 0;
    for(int i=0; i<numThreads; i++){
       average += averageArray[i];
    }

    return average/length;
    
};

double numArrayWrapper::sqrtAverageLogSlow(numArray data){
    double average = 0;
    data.logArray();
    data.sqrtArray();
    for(int i=0; i<data.getLength(); i++){
        average+=data.getElement(i);
    }
    return average/data.getLength(); 
    
};
