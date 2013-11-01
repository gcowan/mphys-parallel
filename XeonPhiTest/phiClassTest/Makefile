CXX		=	icpc
LD		=	icc
DEBUG		=
CXXFLAGS	=	-openmp -O3 -vec-report2
LDFLAGS		=	

EXE		=	phiTest
CXXFILES	=	$(wildcard *.cpp)
OBJECTS		=	$(CXXFILES:.cpp=.o)

$(EXE) : $(OBJECTS)
	$(LD) -V $^ $(CXXFLAGS) $(LDFLAGS) -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :	
	rm -fr $(EXE) $(OBJECTS) $(OBJECTS:.o=MIC.o)
