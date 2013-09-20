#include"numArray.h"
#include <ctime>
#include <iostream>


int main(){
    numArray * a = new numArray(100000000);
    numArray * b = new numArray(100000000);
    //Setting up clock to see how long serial and parallel programms run for
    time_t begining;
    time_t end;
    time(&begining);
    a->serialAdd(*b);
    a->serialAdd(*b);
    a->serialMultiply(2);
    a->serialAdd(*b);
    a->serialMultiply(2);
    a->serialAdd(*b);
    a->serialMultiply(2);
    time(&end);
    std::cout<<"The time that took to run was "<<difftime(end,begining)<<" the last number was "<<a->getElement(99999999)<<std::endl;
    a->deleteNumArray();
    delete a;
    a = new numArray(100000000);
    time(&begining);
    a->parallelAdd(*b);
    a->parallelAdd(*b);
    a->parallelMultiply(2); 
    a->serialAdd(*b);
    a->serialMultiply(2);
    a->serialAdd(*b);
    a->serialMultiply(2);
    time(&end);
    std::cout<<"The time that took to run was "<<difftime(end,begining)<<" the last number was "<<a->getElement(99999999)<<std::endl;
    //have to include these functions first as calling delete[] from the destructor caused unnallocated pointer runtime errors
    a->deleteNumArray();
    b->deleteNumArray();
    delete a,b;
    return 0;



}
