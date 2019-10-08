MACHINE = PRAO
#MACHINE = STAMPEDE2
LDFLAGS =
DEPENDENCIES = Grid.cpp Field.cpp TimeIntegration.cpp 
SOURCES = main.cpp
EXEC    = a.out

##############################################################
## building object list ##

OTHERLIBS=
OTHERINC=

OBJECTS = $(DEPENDENCIES) $(addsuffix .o,$(basename $(SOURCES)))

########################################################################
ifeq ($(MACHINE),PRAO)

COMPILE = mpic++ -g -O0 -std=c++17
CC= mpicc -O0 -std=c++17
CFLAGS = -g -Wall -Wextra

#COMPILE = mpic++ -O0 -std=c++11 -fopenmp
#CC= mpicc -O0 -std=c++11 -fopenmp
#FC= mpif90 -O0 -fopenmp
#CFLAGS = -g -Wall -Wextra 

#FFTW3_DIR=/home/prao/software/install/fftw-3.3.8
#FFTW3_LIB=-L$(FFTW3_DIR)/lib -lfftw3_mpi -lfftw3 -lm
#FFTW3_INC=-I$(FFTW3_DIR)/include

#HYPRE_DIR=/home/prao/software/hypre/src/hypre
#HYPRE_LIB=-L$(HYPRE_DIR)/lib -lHYPRE -lm
#HYPRE_INC=-I$(HYPRE_DIR)/include

endif
########################################################################
ifeq ($(MACHINE),STAMPEDE2)

COMPILE = icpc -O3 -fopenmp
CC= icpc -O0 -fopenmp
CFLAGS = -Wall -Wextra

#COMPILE = icpc  -O0 -fopenmp
#CC= icpc -O0 -fopenmp
#FC= mpif90 -O0 -fopenmp
#CFLAGS = -g -Wall -Wextra

HYPRE_DIR=$(TACC_HYPRE_DIR)
HYPRE_LIB=-L$(TACC_HYPRE_LIB) -lHYPRE -lm
HYPRE_INC=-I$(TACC_HYPRE_INC)

MPI_LIB=-L$(TACC_IMPI_LIB) -lmpi 
MPI_INC=-I$(TACC_IMPI_INC)

endif
########################################################################
LIBS=$(MPI_LIB) $(FFTW3_LIB) $(HYPRE_LIB) 
INCLUDE=$(MPI_INC) $(FFTW3_INC) $(HYPRE_INC)

## To print a variable
##$(info $(CFLAGS))

## build the executable ##
$(EXEC): $(OBJECTS)
	$(COMPILE) $(INCLUDE) $(OBJECTS) -o  $@ $(LDFLAGS) $(LIBS) 

# build cpp object files
%.o: %.cpp
	@echo building $< 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@


# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod *~ \#* 
	@echo 'files cleaned'

