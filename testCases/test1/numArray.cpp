#include<iostream>
#include"numArray.h"
#include <omp.h>
#include <cmath>


numArray::numArray(int i){
    int j;
    array=new double[i];
    length = i;
    for(j=0; j<i; j++){
        array[j]=j;
    }

};

void numArray::deleteNumArray(){
    delete [] array;
};

__attribute__((target(mic))) void numArray::parallelAdd(numArray  otherArray){
    int i;
    //This is the part which parallelizes the for loop
    #pragma omp parallel for
        for(i=0; i<length; i++){
	    	array[i] += otherArray.getElement(i);
	}
     
};


void numArray::serialAdd(numArray  otherArray){
    int i;
        for(i=0; i<length; i++){
	    	array[i] += otherArray.getElement(i);
	}
     
};


void numArray::printArray(){
    int i;
    for(i=0; i<length; i++){
        std::cout<<array[i]<<std::endl;
    }
};

double numArray::getElement(int k){
    if(k>=length){
        std::cout<<"The element being accessed must be less than the number of elements ("<<k<<","<<length<<")"<<std::endl;
	return 0;
    }
    else{
         return array[k];   
 
    }


};

void numArray::serialMultiply(int n){
   int i;
   for(i =0; i<length; i++){
       array[i]*=n;
   }

};

__attribute__((target(mic))) void numArray::parallelMultiply(int n){
    #pragma omp parallel for
        for(int i=0; i<length; i++){
           array[i]*=n;	
	}
};

void numArray::sqrtArray(){
    for(int i=0; i<length; i++){
        array[i] = sqrt(array[i]);
    } 
};


__attribute__((target(mic))) void numArray::sqrtParallelArray(){
 #pragma omp parallel for 
     for(int i=0; i<length; i++){
         array[i] = sqrt(array[i]);
    }
 };

 void numArray::logArray(){
    for(int i=0; i<length; i++){
        array[i] = log(array[i]);
    } 
};

__attribute__((target(mic))) void numArray::logParallelArray(){
 #pragma omp parallel for 
     for(int i=0; i<length; i++){
         array[i] = log(array[i]);
    }
 };

 int numArray::getLength(){
     return length;
 };

