
#ifndef SRC_PHYSICS_LINEARSYSTEM_SYSTEM_DISTRIBUTEDSYSTEM_H_
#define SRC_PHYSICS_LINEARSYSTEM_SYSTEM_DISTRIBUTEDSYSTEM_H_

#include "linearsystem.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.dense.distributed.h"
#include "math/vector.sparse.h"

#include <vector>

namespace espreso {

class DistributedComposer;

struct DistributedAssemblerData: public AssemblerData {
	MatrixCSRDistributed K, M, C, CM;
	VectorsDenseDistributed R, f, x;
	VectorsSparse BC;

	DistributedAssemblerData(): AssemblerData(&K, &M, &C, &CM, &R, &f, &x, &BC) {}

	void print(const Builder *builder, const char* prefix, const char* suffix);
};

struct DistributedSolverData: public SolverData {
	MatrixCSRDistributed K;
	VectorsDenseDistributed R, f, x, y;
	VectorsSparse BC;

	DistributedSolverData(SystemSolver *linearSolver): SolverData(&K, &R, &f, &x, &y, &BC, linearSolver) {}

	void setDirichlet(const Builder *builder);

	void printData(const Builder *builder, const char* prefix);
	void printSolution(const Builder *builder, const char* prefix);

private:
	std::vector<double> _KDirichletValues;
};

}


#endif /* SRC_PHYSICS_LINEARSYSTEM_SYSTEM_DISTRIBUTEDSYSTEM_H_ */
