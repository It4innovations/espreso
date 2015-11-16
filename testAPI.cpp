
#include "libespreso/espreso.h"
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int MPIrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	ESPRESOInit(MPI_COMM_WORLD);

	std::stringstream ss;
	ss << argv[1] << MPIrank;

	esint n;
	esint nelt;
	esint *eltptr;
	esint *eltvar;
	double *values;
	esint nvar;
	esint valueSize = 0;

	std::stringstream KeInfo;
	KeInfo << ss.str() << "/KeInfo.txt";
	std::ifstream file(KeInfo.str().c_str());
	if (file.is_open()) {
		file >> n;
		file >> nelt;
		file >> nvar;
	}
	file.close();

	eltptr = new esint[nelt + 1];
	eltvar = new esint[nvar];

	std::stringstream KeElmPtr;
	KeElmPtr << ss.str() << "/eltptr.txt";
	file.open(KeElmPtr.str().c_str());
	esint index = 0;
	if (file.is_open()) {
		while (file >> eltptr[index++]);
	}
	file.close();

	for (size_t i = 0; i < nelt; i++) {
		valueSize += (eltptr[i + 1] - eltptr[i]) * (eltptr[i + 1] - eltptr[i]);
	}

	std::stringstream KeIndex;
	KeIndex << ss.str() << "/KeIndex.txt";
	file.open(KeIndex.str().c_str());
	index = 0;
	if (file.is_open()) {
		while (file >> eltvar[index++]);
	}
	file.close();

	values = new double[valueSize];

	std::stringstream KeValues;
	KeValues << ss.str() << "/Ke.txt";
	file.open(KeValues.str().c_str());
	index = 0;
	if (file.is_open()) {
		while (file >> values[index++]);
	}
	file.close();

	ESPRESOMatrix K;
	ESPRESOCreateMatrixElemental(n, nelt, eltptr, eltvar, values, &K);

	ESPRESODoubleVector rhs = new ESPRESOStructDoubleVector();
	ESPRESOMap dirichlet = new ESPRESOStructMap();
	ESPRESOIntVector l2g = new ESPRESOStructIntVector();
	ESPRESOIntVector neighbourRanks = new ESPRESOStructIntVector();


	std::stringstream RHS;
	RHS << ss.str() << "/f.txt";
	file.open(RHS.str().c_str());
	index = 0;
	if (file.is_open()) {
		file >> rhs->size;
		rhs->values = new double[rhs->size];
		while (file >> rhs->values[index++]);
	}
	file.close();

	std::stringstream DIRI;
	DIRI << ss.str() << "/BCDIndex.txt";
	file.open(DIRI.str().c_str());
	index = 0;
	if (file.is_open()) {
		file >> dirichlet->size;
		if (dirichlet->size) {
			dirichlet->indices = new esint[dirichlet->size];
			dirichlet->values = new double[dirichlet->size];
			while (file >> dirichlet->indices[index++]);
		}
	}
	file.close();

	std::stringstream DIRV;
	DIRV << ss.str() << "/BCDValue.txt";
	file.open(DIRV.str().c_str());
	index = 0;
	if (file.is_open()) {
		if (dirichlet->size) {
			while (file >> dirichlet->values[index++]);
		}
	}
	file.close();

	std::stringstream L2G;
	L2G << ss.str() << "/l2g.txt";
	file.open(L2G.str().c_str());
	index = 0;
	if (file.is_open()) {
		file >> l2g->size;
		l2g->values = new esint[l2g->size];
		while (file >> l2g->values[index++]);
	}
	file.close();

	std::stringstream RANKS;
	RANKS << ss.str() << "/ranks.txt";
	file.open(RANKS.str().c_str());
	index = 0;
	if (file.is_open()) {
		file >> neighbourRanks->size;
		neighbourRanks->values = new esint[neighbourRanks->size];
		while (file >> neighbourRanks->values[index++]);
	}
	file.close();

	ESPRESOFETIInstance instance;
	//ESPRESOPrepareFETIInstance(&K, &rhs, &dirichlet, &l2g, &neighbourRanks, &instance);

	ESPRESODoubleVector solution = new ESPRESOStructDoubleVector();
	solution->size = rhs->size;
	solution->values = new double[solution->size];
	//ESPRESOSolveFETI(&instance, solution->size, solution->values);

	ESPRESODoubleVector vvv;
	ESPRESOCreateDoubleVector(rhs->size, rhs->values, &vvv);

	ESPRESODestroy(vvv);

	ESPRESOFinalize();

	MPI_Finalize();
}



