ALL:

ifdef TEC360HOME
CFLAGS		 = -I${TEC360HOME}/include/ -DTECIO=1
LIBS		 = ${TEC360HOME}/lib/tecio64.a -lstdc++
else
CFLAGS		 =
LIBS             =
endif
LDFLAGS		 = 
FFLAGS		 =
CPPFLAGS 	 =  
FPPFLAGS         =
LOCDIR		 = 
MANSEC           = SNES

LIBFLAG          =

SOURCEC = bcs.c bmv.c compgeom.c ibm.c ibm_io.c init.c \
          main.c metrics.c poisson.c rhs.c rheology.c\
	  variables.c fsi.c implicitsolver.c\
	  fsi_move.c solvers.c copepod.c fish.c cstart.c spline.c\
          les.c k-omega.c wallfunction.c rhs2.c poisson_hypre.c platlet.c

OBJSC =  bcs.o bmv.o compgeom.o ibm.o ibm_io.o init.o \
         main.o metrics.o poisson.o rhs.o rheology.o\
         variables.o fsi.o implicitsolver.o\
         fsi_move.o solvers.o copepod.o fish.o cstart.o spline.o\
         les.o k-omega.o wallfunction.o rhs2.o poisson_hypre.o platlet.o\

LIBBASE = libpetscmat

#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/variables
#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/rules

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

testt: ${OBJSC}
	-$(CLINKER) -o testt ${OBJSC} ${PETSC_LIB}

	rm main.o

data: variables.o  compgeom.o data_ibm.o ibm_io.o fsi.o fsi_move.o fish.o  data.o
	-${CLINKER} -o data variables.o compgeom.o data_ibm.o ibm_io.o fsi.o fsi_move.o fish.o data.o ${PETSC_LIB} ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}

itfcsearch: itfcsearch.o variables.o compgeom.o
	-${CLINKER} -o itfcsearch itfcsearch.o  variables.o compgeom.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

data_vtk: data_surface.o 
	-${CLINKER} -o data_vtk data_surface.o  ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

datalis: data_file2lis.o
	-${CLINKER} -o datalis data_file2lis.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}

datafile: data_list2file.o
	-${CLINKER} -o datafile data_list2file.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}
cleanobj:
	rm -f *.o

include ${PETSC_DIR}/lib/petsc/conf/test
#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/test
