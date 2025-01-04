F77 = gfortran


OBJ = JetOutput-0.8.o Calculate.o
all: Calculate.exe

JetOutput-0.8.o: JetOutput-0.8.for
	$(F77) -c -fno-align-commons JetOutput-0.8.for

Calculate.o: Calculate.for
	$(F77) -c Calculate.for

Calculate.exe: $(OBJ)	
	$(F77) -O $^ -o $@


clean:
	rm -f *.o *.*~ *.exe
