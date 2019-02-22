CXX = g++
CXXFLAGS = -Wall -fopenmp -O3
PROGS = wvstrain

.PHONY: all clean

all: $(PROGS)

clean:
	rm -f $(PROGS) *.o

wvstrain: krakendb.o quickfile.o krakenutil.o 

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp
