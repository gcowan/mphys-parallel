CC=icpc
LDFLAGS= -mkl 

CFLAGS=-openmp -openmp-report2 -vec-report3 

ROOTFLAGS=-L/data/s1007232/root/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread  -ldl -rdynamic -pthread -m64 -I/data/s1007232/root/root/include -lMinuit

prog : main.o myGauss.o 
	$(CC) $(ROOTFLAGS) $(CFLAGS) $(LDFLAGS)  main.o myGauss.o -o prog -lm -lpthread

main.o : main.cpp
	$(CC)  -c $(ROOTFLAGS)  -vec-report2 -std=c++11 -openmp-report2 main.cpp -lm -lpthread

myGauss.o : myGauss.cpp
	$(CC) $(CFLAGS) -c   myGauss.cpp -lm -lpthread

clean :
	rm prog *.o
