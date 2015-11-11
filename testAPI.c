
#include "libespreso/espreso.h"

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	ESPRESOInit(MPI_COMM_WORLD);

	ESPRESOMat stiffnessMatrix;

	ESPRESODoubleVector rhs;
	rhs.size = 5;
	rhs.values = (double*)malloc(5 * sizeof(double));

	ESPRESOMap dirichlet;
	dirichlet.size = 4;
	dirichlet.indices = (int*)malloc(5 * sizeof(double));
	dirichlet.values = (double*)malloc(5 * sizeof(double));

	ESPRESOIntVector l2g;
	l2g.size = 5;
	l2g.values = (int*)malloc(5 * sizeof(double));

	ESPRESOIntVector neighbourRanks;
	neighbourRanks.size = 3;
	neighbourRanks.values = (int*)malloc(3 * sizeof(double));

	ESPRESODoubleVector solution;
	solution.size = 5;
	solution.values = (double*)malloc(5 * sizeof(double));

	ESPRESOSolveFETI(&stiffnessMatrix, &rhs, &dirichlet, &l2g, &neighbourRanks, &solution);

	ESPRESOFree(stiffnessMatrix);

	MPI_Finalize();
}



