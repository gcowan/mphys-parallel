CC=icpc

CFLAGS=-openmp -openmp-report2 -vec-report2 -O3

FITFLAGS=  -openmp -mkl 



LINK =   -mkl -lpthread -lm 


prog : mainProg.cpp  
	$(CC)  $(CFLAGS) $(FITFLAGS)  -c fit.cpp 
	$(CC) $(CFLAGS) -c mainProg.cpp
	$(CC)  -o prog  mainProg.o fit.o  $(LINK)

clean :
	rm prog *.o
