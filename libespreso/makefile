

all:
	mpiicpc -std=c++11 apiexample.cpp -L. -lfeti4i -L../libs -lmetis32 -o apiexample

static:
	mpiicpc -mkl -std=c++11 apiexample.cpp -Wl,-Bstatic -L. -lfeti4i -L../libs -lmetis32 -Wl,-Bdynamic
