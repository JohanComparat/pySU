#FC = gfortran 
#FFLAGS =   -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
#LDFLAGS =  -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
FC = ifort 
FFLAGS =   -O2  -g -traceback -ftz -unroll  -openmp  -i-dynamic -mcmodel=medium -convert big_endian
LDFLAGS =  -O2  -g -traceback -ftz -unroll  -openmp  -i-dynamic -mcmodel=medium -convert big_endian
#FC     = pgf77
#FFLAGS = -O3  -mp  -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee
#LDFLAGS = -O3  -mp   -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee

OBJ = PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2mod_analyze.o

PMP2main: $(OBJ)  PMP2main.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2analyze: $(OBJ)  PMP2analyze.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2init: PMP2init.o
	$(FC) $(LDFLAGS)  -o $@.exe $^
PMP2start: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2start.o
	$(FC) $(LDFLAGS) -o $@.exe $^
.f.o: 
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(FC) $(FFLAGS) -w -c $<

clean:
	rm -f *.o *.exe
