CXX     = g++
MPICXX  = mpic++
LIB     = -lmpi -lpthread -lgfortran -llapack -lblas
#LIB    += -L./lib 
FLAG    = -Wall -I$(INCLUDE) -O2
DEF_OPENMP = -fopenmp -DCAPLET_OPENMP
DEF_MPI = -DCAPLET_MPI

SRC        = src
OBJ_MPI    = obj/mpi
OBJ_OPENMP = obj/openmp
BIN        = bin
TMP        = tmp
INCLUDE    = include
EXAMPLE    = example

CAPLET_MPI_OBJ = \
	$(OBJ_MPI)/caplet.o \
	$(OBJ_MPI)/caplet_elem.o \
	$(OBJ_MPI)/caplet_int.o \
	$(OBJ_MPI)/caplet_widgets.o \
	$(OBJ_MPI)/main.o \

CAPLET_OPENMP_OBJ = \
	$(OBJ_OPENMP)/caplet.o \
	$(OBJ_OPENMP)/caplet_elem.o \
	$(OBJ_OPENMP)/caplet_int.o \
	$(OBJ_OPENMP)/caplet_widgets.o \
	$(OBJ_OPENMP)/main.o \

all: capletMPI capletOpenMP copy

mpi: capletMPI copy

openmp: capletOpenMP copy

capletMPI: $(CAPLET_MPI_OBJ)
	$(MPICXX) $(LIB) $(FLAG) $(DEF_MPI) -o $(BIN)/$@ $^ 

$(OBJ_MPI)/%.o: $(SRC)/%.cpp
	$(MPICXX) $(LIB) $(FLAG) $(DEF_MPI) -c $< -o $@

capletOpenMP: $(CAPLET_OPENMP_OBJ)
	$(MPICXX) $(LIB) $(FLAG) $(DEF_OPENMP) -o $(BIN)/$@ $^

$(OBJ_OPENMP)/%.o: $(SRC)/%.cpp
	$(MPICXX) $(LIB) $(FLAG) $(DEF_OPENMP) -c $< -o $@

copy: 
	cp $(BIN)/* $(TMP); cp $(EXAMPLE)/* $(TMP)

.phony: clean
clean:
	rm -rf $(OBJ_MPI)/* $(OBJ_OPENMP)/* $(TMP)/*