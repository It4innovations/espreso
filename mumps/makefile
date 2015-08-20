.DEFAULT_GOAL=all

BINARIES=mumps_solve
CC=gcc
FC=gfortran
STD=-std=c11

HOME_DIR=/home/bes0030

MUMPS_DIR=${HOME_DIR}/MUMPS_5.0.1
MUMPS_LIB=${MUMPS_DIR}/lib
MUMPS_LIBSEQ=${MUMPS_DIR}/libseq
MUMPS_INC=${MUMPS_DIR}/include

PORD_DIR=${MUMPS_DIR}/PORD
PORD_LIB=${PORD_DIR}/lib

METIS_DIR=${HOME_DIR}/metis-5.1.0
METIS_LIB=${METIS_DIR}/build/Linux-x86_64/libmetis

SCOTCH_DIR=${HOME_DIR}/scotch_6.0.4
SCOTCH_LIB=${SCOTCH_DIR}/lib
SCOTCH_FLAGS=-lesmumps -lscotch -lscotcherr

BLAS_LIB=${MUMPS_DIR}
BLAS_FLAGS=-lblas

.PHONY=all
all: mumps_solve

mumps_solve: mumps_solve.o
	${FC} mumps_solve.o -o mumps_solve -L${MUMPS_LIB} -ldmumps ${MUMPS_LIB}/libmumps_common.a ${MUMPS_LIBSEQ}/libmpiseq.a ${METIS_LIB}/libmetis.a -L${PORD_LIB} -lpord -L${SCOTCH_LIB} ${SCOTCH_FLAGS} -L${BLAS_LIB} ${BLAS_FLAGS} -lpthread 

mumps_solve.o: mumps_solve.c
	${CC} -c mumps_solve.c -o mumps_solve.o -I${MUMPS_LIBSEQ} -I${MUMPS_INC} ${STD}

.PHONY=clean
clean:
	rm -f ${BINARIES} *.o
