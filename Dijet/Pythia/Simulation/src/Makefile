CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IPYTHIA = ${shell pythia8-config --cxxflags}
LPYTHIA = ${shell pythia8-config --libs}
# IPYTHIA = -I/Users/duanpanpan/Desktop/Work/Model/PYTHIA8/pythia8310/include
# LPYTHIA = -L/Users/duanpanpan/Desktop/Work/Model/PYTHIA8/pythia8310/lib -lpythia8 -ldl


TARGET = Calculate.exe
all: $(TARGET)

Pythia_Multiplicity_ATLAS.o: Pythia_Multiplicity_ATLAS.cpp
	$(CXX) $(CXXFLAGS) $(IPYTHIA) -c Pythia_Multiplicity_ATLAS.cpp

$(TARGET): Pythia_Multiplicity_ATLAS.o
	$(CXX) $(CXXFLAGS) $^ $(LPYTHIA) -o $@


clean:
	rm -f *.o *.*~ *.exe
