FC := mpif90
CC := mpicc

build: test

feti4i_mod.o: ../feti4i_mod.F90
	${FC} -c $^

feti4i.mod: 	../feti4i_mod.F90 feti4i_mod.o
	@true

test.o: test.F90 feti4i.mod
	${FC} -c -g test.F90 -o $@

test: test.o
	${LINK.F} test.o -o $@ -L.. -lfeti4i

clean:
	${RM} test.o feti4i_mod.o feti4i.mod

.PHONY: build clean

