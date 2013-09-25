#include <iostream>
#include "numArray.h"
#include "numArrayWrapper.h"
#include <chrono>
#include <fstream>

int main(){
    std::ofstream output;
    output.open("arrayValues.txt");
    int arraySize = 1;
    numArray * a = new numArray[arraySize];
    numArrayWrapper wrapper;
    //Initialising the array
    for(int i=0; i<arraySize; i++){
        a[i].generateNumArray(10000);
    }
    //creating the timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //now doing slow sqrt
    double totalAverage = 0;
    output<<"Slow averages"<<std::endl; 
    for(int i=0; i<arraySize; i++){
        double val = wrapper.sqrtAverageLogSlow(a[i]);
	totalAverage+=val;
	output<<val<<std::endl;
	
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end-start;
    std::cout<<"The value of the totalAverage done by serial computation is: "<<totalAverage/arraySize<<" the calculation took "<<elapsedSeconds.count()<<std::endl;

    //Now cleaning up array
    for(int i=0; i<arraySize; i++){
        a[i].deleteNumArray();
    }

    delete[] a;

    //reinitializing array
    a = new numArray[arraySize];
    for(int i=0; i<arraySize; i++){
        a[i].generateNumArray(10000);
    }

    start = std::chrono::system_clock::now();
    output<<"Fast Averages"<<std::endl;
    //now doing fast sqrt
    totalAverage = 0;
    for(int i=0; i<arraySize; i++){
        double val = wrapper.sqrtAverageLogFast(a[i]);
	totalAverage += val;
	output<<val<<std::endl;
    }
    end = std::chrono::system_clock::now();
    elapsedSeconds = end-start;
    std::cout<<"The value of the totalAverage done by parallel computation is: "<<totalAverage/arraySize<<" the calculation took "<<elapsedSeconds.count()<<std::endl;
    return 0; 


    

}
