dir=class
OPT=-O3 -zero -qopenmp -mkl=parallel  
FFLAGS= $(OPT) -L$(dir)/lib/ -I$(dir)/mod/  
F90 = ifort
all: flow.out clean
libg:
	cd class/;make;make install
libs:
	rm -f *.out
	cd class/; make;make install
%.out:%.f90 
	$(F90) $(FFLAGS) $<  -lgologie  -o  $@
%.o: %.f90 
	$(F90) $(FFLAGS) -c $<
clean:
	rm -f  *.o *.mod *~
cleanout:
	rm -f *.out
cleantout: clean cleanout
	cd class/; make clean
