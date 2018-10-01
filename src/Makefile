CXX = g++
CXXFLAGS = -Wall -fopenmp -O3
PROGS = mladen_main 

.PHONY: all clean

all: $(PROGS)

clean:
	rm -f $(PROGS) *.o

mladen_main: krakendb.o quickfile.o krakenutil.o 

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp
