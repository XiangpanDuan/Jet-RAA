all:
	gfortran -c -fno-strict-overflow -fno-align-commons JetOutput-0.8.for vint.f90 vish_gen.f90 main.f90
	gfortran -o generate.exe *.o


clean:
	rm -f *.o *.*~ generate.exe *.mod