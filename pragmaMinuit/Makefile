CC=icpc

CFLAGS=-openmp -openmp-report2 -vec-report2 -O3

FITFLAGS=  -openmp -mkl 

ROOTFLAGS =  -pthread -m64 -I/data/s1007232/root/root/include 


LINK =   -mkl -L/data/s1007232/root/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -lMinuit -rdynamic  -lpthread -lm 


prog : mainProg.cpp  
	$(CC)  $(CFLAGS) $(FITFLAGS)  -c fit.cpp 
	$(CC) $(CFLAGS) $(ROOTFLAGS) -c mainProg.cpp
	$(CC)  -o prog  mainProg.o fit.o  $(LINK)

clean :
	rm prog *.o
