CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

# IGSL = ${shell gsl-config --cflags}
# LGSL = ${shell gsl-config --libs}
IGSL = -I/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/include
LGSL = -L/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = CalOmegaC.exe
all: $(TARGET)

TrilinearInterpolation.o: TrilinearInterpolation.cpp TrilinearInterpolation.h
	$(CXX) $(CXXFLAGS) $(IGSL) -c TrilinearInterpolation.cpp

CalOmegaC.o: CalOmegaC.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) -c CalOmegaC.cpp

$(TARGET): CalOmegaC.o TrilinearInterpolation.o
	$(CXX) $(CXXFLAGS) $^ $(LGSL) -o $@


clean:
	rm -f *.o *.*~ *.exe
