#include"numArray.h"
#include <ctime>
#include <iostream>


int main(){
    numArray * a = new numArray(100000);
    numArray * b = new numArray(100000);
    a->parallelAdd(*b);
    a->printArray();
    //delete a,b;
    return 0;



}
