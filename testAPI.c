
#include "libespreso/espreso.h"

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	ESPRESO_Init(MPI_COMM_WORLD);
	ESPRESO_Finalize();

	LocalStiffnessMatrices Ke;
	Ke.size = 5;
	Ke.array = (int*)malloc(sizeof(int) * 5);
	int i;
	for (i = 0; i < 5; i++) {
		Ke.array[i] = i;
	}

	ESPRESOHandler elasticity;
	DoubleVector rhs;
	DoubleVector solution;
	FETI_PrepareElasticity(&Ke, &rhs, &elasticity);
	Solve(&elasticity, &solution);

	MPI_Finalize();
}



