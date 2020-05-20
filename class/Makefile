
#F90 = /home/asedeki/intel/fc/9.1.043/bin/ifort

F90 = ifort
OBJ = parameters.o class_system.o ode.o rkf45.o class_indice.o  class_bublesC.o class_interaction.o class_responsefunction.o class_FreeEnergie.o class_quasi1D.o
mod = *.mod
dir=.
OPT=-O3 -zero -qopenmp -mkl=parallel  
FFLAGS2=$(OPT) -L$(dir)/lib/ -I$(dir)/mod
all: lib

lib: $(OBJ)	
	ar r libgologie.a $(OBJ) ; ranlib libgologie.a

%.o: %.f90 
	$(F90) $(FFLAGS2) -c $<
install:
	mv libgologie.a $(dir)/lib/
	cp $(mod)  $(dir)/mod
clean: 
	rm -f *.mod *.o 

