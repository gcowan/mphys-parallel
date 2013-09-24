#include"numArray.h"
#include <ctime>
#include <iostream>


int main(){
    numArray * a = new numArray(100000000);
    //Setting up clock to see how long serial and parallel programms run for
    time_t begining;
    time_t end;
    time(&begining);
    a->sqrtArray();
    a->logArray();
    time(&end);
    std::cout<<"The time that took to run was "<<difftime(end,begining)<<" the last number was "<<a->getElement(99999999)<<std::endl;
    a->deleteNumArray();
    delete a;
    a = new numArray(100000000);
    int length = a->getLength();
    std::cout<<length<<std::endl;
    time(&begining);
    #pragma offload target(mic) in(b:length(length)) inout(a:length(length))
    a->sqrtParallelArray();
    a->logParallelArray();
    time(&end);
    std::cout<<"The time that took to run was "<<difftime(end,begining)<<" the last number was "<<a->getElement(99999999)<<std::endl;
    //have to include these functions first as calling delete[] from the destructor caused unnallocated pointer runtime errors
    a->deleteNumArray();
    
    delete a;
    return 0;



}
