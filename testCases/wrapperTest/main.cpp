#include <iostream>
#include "numArray.h"
#include "numArrayWrapper.h"
#include <chrono>

int main(){
    int arraySize = 25;
    numArray * a = new numArray[arraySize];
    numArrayWrapper wrapper;
    //Initialising the array
    for(int i=0; i<arraySize; i++){
        a[i].generateNumArray(100000);
    }
    //creating the timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //now doing slow sqrt
    double totalAverage = 0;
    for(int i=0; i<arraySize; i++){
        totalAverage+= wrapper.sqrtAverageSlow(a[i]);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end-start;
    std::cout<<"The value of the totalAverage done by serial computation is: "<<totalAverage/arraySize<<" the calculation took "<<elapsedSeconds.count()<<std::endl;

    //Now cleaning up and reinitializing array
    for(int i=0; i<arraySize; i++){
        a[i].deleteNumArray();
	a[i].generateNumArray(100000);
    }

    start = std::chrono::system_clock::now();
    //now doing slow sqrt
    totalAverage = 0;
    for(int i=0; i<arraySize; i++){
        totalAverage+= wrapper.sqrtAverageFast(a[i]);
    }
    end = std::chrono::system_clock::now();
    elapsedSeconds = end-start;
    std::cout<<"The value of the totalAverage done by parallel computation is: "<<totalAverage/arraySize<<" the calculation took "<<elapsedSeconds.count()<<std::endl;
    return 0; 


    

}
