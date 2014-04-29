PLAT   = _pdsyev_darwin
FOR    = mpif90
CC = icc

FFLAGS = -O3 -ipo -xHost -traceback -g
#-align -mcmodel=medium -traceback

SCALAPACK = -lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 -mkl=sequential -lmpi

#LIB    =  $(SCALAPACK)
#LIB    =  -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -lpthread -lm

LIB = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#LIB = -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_sequential -lmkl_core -lpthread

OBJ = accuracy.o timer.o input.o pdlaprnt_local.o  
#utils.o

default:   diag$(PLAT).x

diag$(PLAT).x:  $(OBJ) diag_pdsyev.o
	$(FOR) $(FFLAGS) -o diag$(PLAT).x $(OBJ) diag_pdsyev.o $(LIB)

diag_pdsyev.o:  diag_pdsyev.f90 $(OBJ) 
	$(FOR) -c diag_pdsyev.f90 $(FFLAGS)

pdlaprnt_local.o: pdlaprnt_local.f
	$(FOR) -c pdlaprnt_local.f $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90 accuracy.o
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90 accuracy.o
	$(FOR) -c input.f90 $(FFLAGS)

#.c.o:
#	$(CC) -openmp -c -O2 -g -DUSE_COMPLEX $*.c


clean:
	rm $(OBJ) *.mod diag_pdsyev.o *.x

