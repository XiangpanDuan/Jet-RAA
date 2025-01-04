CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IROOT = ${shell root-config --cflags}
LROOT = ${shell root-config --libs --glibs}
IFASTJET = ${shell fastjet-config --cxxflags --plugins}
LFASTJET = ${shell fastjet-config --libs --plugins}
# IFASTJET = ${shell /Users/duanpanpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --plugins}
# LFASTJET = ${shell /Users/duanpanpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --libs --plugins}
# FASTJET = `/Users/duanpanpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
# IFASTJET = -I/Users/duanpanpan/Desktop/Work/Model/Fastjet/fastjet-install/include
# LFASTJET = -L/Users/duanpanpan/Desktop/Work/Model/Fastjet/fastjet-install/lib -lfastjettools -lfastjet


TARGET = Analysis.exe
all: $(TARGET)

Jet_Multiplicity_Dijet_ATLAS_8000GeV.o: Jet_Multiplicity_Dijet_ATLAS_8000GeV.cpp
	$(CXX) $(CXXFLAGS) $(IROOT) $(IFASTJET) -c Jet_Multiplicity_Dijet_ATLAS_8000GeV.cpp

$(TARGET): Jet_Multiplicity_Dijet_ATLAS_8000GeV.o
	$(CXX) $(CXXFLAGS) $^ $(LROOT) $(LFASTJET) -o $@


clean:
	rm -f *.o *.*~ *.exe
