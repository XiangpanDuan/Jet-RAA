F77 = gfortran


TARGET = Analysis.exe
all: $(TARGET)

JetOutput-0.8.o: JetOutput-0.8.for
	$(F77) -c -fno-align-commons JetOutput-0.8.for

Analysis.o: Analysis.for
	$(F77) -c Analysis.for

$(TARGET): Analysis.o JetOutput-0.8.o
	$(F77) -O $^ -o $@


clean:
	rm -f *.o *.*~ *.exe
