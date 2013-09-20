#ifndef NUMARRAY_H
#define NUMARRAY_H
//Simple class to store an array of numbers
class numArray{
   private: 
         int * array;
         int length;

    public:
    	 //Constructor constructs ints of length n
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
	 int getElement(int i);
	 //Multiply every element by a number
	 void serialMultiply(int n);
	 void parallelMultiply(int n);
};

#endif
