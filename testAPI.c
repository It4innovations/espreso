
#include "libespreso/espreso.h"

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	ESPRESOInit(MPI_COMM_WORLD);

	esint n = 5;
	esint nelt = 2;
	esint nvar = 6;

	esint *eltptr = (esint*)malloc(sizeof(esint) * 3);
	eltptr[0] = 1;
	eltptr[1] = 4;
	eltptr[2] = 7;

	esint *eltvar = (esint*)malloc(sizeof(esint) * 6);
	eltvar[0] = 1;
	eltvar[1] = 2;
	eltvar[2] = 3;
	eltvar[3] = 3;
	eltvar[4] = 4;
	eltvar[5] = 5;

	double *values = (double*)malloc(sizeof(double) * 18);
	values[0] = -1;
	values[1] = 2;
	values[2] = 1;
	values[3] = 2;
	values[4] = 1;
	values[5] = 1;
	values[6] = 3;
	values[7] = 1;
	values[8] = 1;
	values[9] = 2;
	values[10] = 1;
	values[11] = 3;
	values[12] = -1;
	values[13] = 2;
	values[14] = 2;
	values[15] = 3;
	values[16] = -1;
	values[17] = 1;

	ESPRESOMat K;

	ESPRESOCreateMatrixElemental(n, nelt, eltptr, eltvar, values, &K);

//	ESPRESOMat stiffnessMatrix;
//
//	ESPRESODoubleVector rhs;
//	rhs.size = 5;
//	rhs.values = (double*)malloc(5 * sizeof(double));
//
//	ESPRESOMap dirichlet;
//	dirichlet.size = 4;
//	dirichlet.indices = (int*)malloc(5 * sizeof(double));
//	dirichlet.values = (double*)malloc(5 * sizeof(double));
//
//	ESPRESOIntVector l2g;
//	l2g.size = 5;
//	l2g.values = (int*)malloc(5 * sizeof(double));
//
//	ESPRESOIntVector neighbourRanks;
//	neighbourRanks.size = 3;
//	neighbourRanks.values = (int*)malloc(3 * sizeof(double));
//
//	ESPRESODoubleVector solution;
//	solution.size = 5;
//	solution.values = (double*)malloc(5 * sizeof(double));
//
//	ESPRESOSolveFETI(&stiffnessMatrix, &rhs, &dirichlet, &l2g, &neighbourRanks, &solution);
//
//	ESPRESOFree(stiffnessMatrix);

	ESPRESOFinalize();

	MPI_Finalize();
}



