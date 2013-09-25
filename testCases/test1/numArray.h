#ifndef NUMARRAY_H
#define NUMARRAY_H
//Simple class to store an array of numbers
class numArray{
   private: 
         double * array;
         int length;

    public:
    	 //Constructor constructs array of doubles of length n
         numArray(int i);
	 //Destructor to remove the array
	 void deleteNumArray();
	 //Method to parallel add
	 void parallelAdd(numArray otherArray);
	 //Method to serial add
	 void serialAdd(numArray otherArray);
	 //Method to print array to screen
	 void printArray();
	 //Method to get ith element of array
	 double getElement(int i);
	 //Multiply every element by a number
	 void serialMultiply(int n);
	 void parallelMultiply(int n);
	 //Method to squareRoot every element in the array
	 void sqrtArray();
	 void sqrtParallelArray();
	 //Method to log every element in the array
	 void logArray();
	 void logParallelArray();
	 //method to get length of array
	 int getLength();
};

#endif
