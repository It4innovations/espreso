CC=mpicc

USRCFLAGS=
DEBUGCFLAGS=-O0 -g -cilk-serialize 

CFLAGS =-c -fopenmp -mkl=parallel -I. -I../solver/ -I./include -DHFETI_SOLVER #-Wall
opt   : CFLAGS += -O3 -xHost -ip -g
debug : CFLAGS += $(DEBUGCFLAGS)
pipe  : CFLAGS += $(USRCFLAGS)
mic   : CFLAGS += -mmic 

LDFLAGS = -mkl=parallel -L. -lsolver -ltbb
LDFLAGS_MIC =  -mmic -mkl=parallel -L. -lsolver_mic -ltbb

SOURCES = Assembler.cpp BoundaryCondition.cpp BoundaryCondOnFace.cpp BSparse.cpp BSparseC.cpp Cluster_g.cpp ContactBC.cpp Coordinate.cpp CoupleIntDouble.cpp CustomData.cpp Data.cpp DirichletBC.cpp Domain_g.cpp ElementLibrary.cpp Fem.cpp KSparse.cpp LinearAlgebra.cpp main.cpp Mesh.cpp NeumannBC.cpp RectangularFace.cpp RSparse.cpp SaveFEM.cpp Solid45.cpp Solid45EqNumbGlobal.cpp Solid45NodNumbLoc.cpp Solver.cpp StiffnessLocal.cpp Test.cpp utility.cpp

OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLE = pmcube
DEBUG_EXEC = pmcube_db
MIC_EXEC   = pmcube_mic
PIPE_EXEC  = pmcube_pipe

opt: $(SOURCES) $(EXECUTABLE)

pipe: $(SOURCES) $(PIPE_EXEC)

debug: $(SOURCES) $(DEBUG_EXEC)

mic: $(SOURCES) $(MIC_EXEC)


$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(DEBUG_EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(MIC_EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS_MIC) $(OBJECTS) -o $@

$(PIPE_EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o 

cleanall:
	rm -rf *.o $(EXECUTABLE) $(DEBUG_EXEC) $(PIPE_EXEC) $(MIC_EXEC)
