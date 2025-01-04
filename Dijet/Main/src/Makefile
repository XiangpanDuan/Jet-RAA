CXX = /opt/homebrew/Cellar/llvm/18.1.8/bin/clang++  #reset g++ for OMP
CXXFLAGS = -std=c++11 -Wall -O3 -fopenmp  #OMP:Open Multi-Processing

IOMP = -I/opt/homebrew/Cellar/llvm/18.1.8/include
LOMP = -L/opt/homebrew/Cellar/llvm/18.1.8/lib
IPDF = -I/opt/homebrew/Cellar/lhapdf/6.5.3/include
LPDF = -L/opt/homebrew/Cellar/lhapdf/6.5.3/lib -lLHAPDF
IGSL = -I/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/include
LGSL = -L/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm
#IPDF = ${shell lhapdf-config --cflags}
#LPDF = ${shell lhapdf-config --libs}
#IGSL = ${shell gsl-config --cflags}
#LGSL = ${shell gsl-config --libs}


TARGET = AA.exe
all: $(TARGET)


Particle.o: Particle.cpp Particle.h
	$(CXX) $(CXXFLAGS) -c Particle.cpp

# ElectroWeak.o: ElectroWeak.cpp ElectroWeak.h
# 	$(CXX) $(CXXFLAGS) $(IGSL) -c ElectroWeak.cpp

QCD.o: QCD.cpp QCD.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c QCD.cpp

DiJetLO.o: DiJetLO.cpp DiJetLO.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c DiJetLO.cpp

EnergyLossOmegaC.o: EnergyLossOmegaC.cpp EnergyLossOmegaC.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c EnergyLossOmegaC.cpp

EnergyLoss.o: EnergyLoss.cpp EnergyLoss.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c EnergyLoss.cpp

Glauber.o: Glauber.cpp Glauber.h
	$(CXX) $(CXXFLAGS) $(IGSL) -c Glauber.cpp

DiJetLOELoss.o: DiJetLOELoss.cpp DiJetLOELoss.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c DiJetLOELoss.cpp


#Calculate cross section
#pp
CalDiJetLO.o: CalDiJetLO.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c CalDiJetLO.cpp

pp.exe: CalDiJetLO.o DiJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@

#PbPb
CalDiJetLOELoss.o: CalDiJetLOELoss.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c CalDiJetLOELoss.cpp

AA.exe: CalDiJetLOELoss.o DiJetLOELoss.o Glauber.o EnergyLoss.o EnergyLossOmegaC.o DiJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@

#PbPb Fixed
CalDiJetLOFixedELoss.o: CalDiJetLOFixedELoss.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c CalDiJetLOFixedELoss.cpp

AAFixed.exe: CalDiJetLOFixedELoss.o DiJetLOELoss.o Glauber.o EnergyLoss.o EnergyLossOmegaC.o DiJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@

#Test
Test.o: Test.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c Test.cpp

Test.exe: Test.o DiJetLOELoss.o Glauber.o EnergyLoss.o EnergyLossOmegaC.o DiJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@

# #All
# CalDiJetLOELoss.o: CalDiJetLOELoss.cpp
# 	$(CXX) $(CXXFLAGS) $(IOMP) $(IPDF) $(IGSL) -c CalDiJetLOELoss.cpp

# $(TARGET): CalDiJetLOELoss.o DiJetLOELoss.o Glauber.o EnergyLoss.o EnergyLossOmegaC.o DiJetLO.o QCD.o Particle.o
# 	$(CXX) $(CXXFLAGS) $^ $(LOMP) $(LPDF) $(LGSL) -o $@



clean:
	rm -f *.o *.*~ *.exe
