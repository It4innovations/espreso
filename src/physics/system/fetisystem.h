
#ifndef SRC_PHYSICS_LINEARSYSTEM_SYSTEM_FETISYSTEM_H_
#define SRC_PHYSICS_LINEARSYSTEM_SYSTEM_FETISYSTEM_H_

#include "linearsystem.h"
#include "feti/generic/FETISystemSolver.h"
#include "math/matrix.csr.feti.h"
#include "math/matrix.ijv.feti.h"
#include "math/matrix.dense.feti.h"
#include "math/vector.dense.feti.h"
#include "math/vector.sparse.h"

#include <vector>

namespace espreso {

struct FETIAssemblerData: public AssemblerData {
	MatrixCSRFETI K, M, C, CM;
	VectorsDenseFETI R, f, x;
	VectorsSparse BC;
	VectorsSparse gapDirection, gap;
	MatrixIJVFETI mortars;

	MatrixIJVFETI B0;
	MatrixCSRFETI RegMat;
	MatrixDenseFETI N1, N2;

	FETIAssemblerData(): AssemblerData(&K, &M, &C, &CM, &R, &f, &x, &BC) {}

	void print(const Builder *builder, const char* prefix, const char* suffix);
};

struct FETISolverData: public SolverData {
	MatrixCSRFETI K;
	VectorsDenseFETI R, f, x, y;
	VectorsSparse BC;
	VectorsSparse gapDirection, gap;
	MatrixIJVFETI mortars;

	MatrixCSRFETI origK, RegMat;
	MatrixDenseFETI N1, N2;
	MatrixIJVFETI B1Dirichlet, B1Gluing, B1Inequality, B0;
	VectorDenseFETI B1c, B1duplication, B1gap;
	VectorDenseFETI Kdiag;
	std::vector<esint> B1Map;

	FETISystemSolver solver;

	FETISolverData(FETIConfiguration &configuration)
	: SolverData(&K, &R, &f, &x, &y, &BC, &solver),
	  solver(configuration, *this) { }

	void buildB1();
	void buildB0();
	void setDirichlet(const Builder *builder);

	bool hasReactionForces() { return true; }
	void printData(const Builder *builder, const char* prefix);
	void printSolution(const Builder *builder, const char* prefix);
};

struct FETISystem: public LinearSystem {

	virtual int nassemblers() { return assemblers.size(); }
	virtual int nsolvers() { return solvers.size(); }

	std::vector<FETIAssemblerData> assemblers;
	std::vector<FETISolverData> solvers;

	FETIAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
	FETISolverData* solver(int index = 0) { return solvers.data() + index; }

	FETISystem(int assemblers, int solvers, FETIConfiguration &configuration);

protected:
	void _builderInit();
	void _builderReset();
	void _builderCreateSystem();
	void _builderUpdateSolution();
};

}



#endif /* SRC_PHYSICS_LINEARSYSTEM_SYSTEM_FETISYSTEM_H_ */
