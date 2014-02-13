CC=icpc
LDFLAGS= -mkl -offload-option,mic,compiler,"-z defs"

CFLAGS=-openmp -openmp-report2 -vec-report2 -g

ROOTFLAGS=-L/data/s1007232/root/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -pthread -m64 -I/data/s1007232/root/root/include -lMinuit

prog : main.o myGauss.o 
	$(CC) $(ROOTFLAGS) $(CFLAGS) $(LDFLAGS) main.o myGauss.o -o prog -lm -lpthread

main.o : main.cpp
	$(CC) -g -c $(ROOTFLAGS) -vec-report2 -openmp-report2 main.cpp

myGauss.o : myGauss.cpp
	$(CC) $(CFLAGS) -c   myGauss.cpp -lm -lpthread

clean :
	rm prog *.o
