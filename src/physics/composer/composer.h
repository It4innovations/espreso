
#ifndef SRC_PHYSICS_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_COMPOSER_H_

#include "physics/system/builder/builder.h"
#include <cstddef>
#include <vector>

namespace espreso {

class Kernel;
class KernelOpt;
struct AssemblerData;
struct IJ;
class Vectors;
class VectorSparse;
class SolverDataProvider;

class Composer {

public:
	Composer(Kernel *kernel, KernelOpt *opt);

	virtual void init() = 0;
	virtual void assemble(const Builder &builder) = 0;

	virtual int esize(esint interval) =0;
	virtual int bsize(esint region, esint interval) =0;
	SolverDataProvider* provider();
	int solutions();

	void solutionChanged(Vectors *solution);
	void nextSubstep();
	void processSolution();

	Kernel *kernel;
	KernelOpt *opt;

	virtual ~Composer();
protected:
	esint getMatrixSize(esint size, bool omitLower);
	void insertKPattern(IJ *target, const esint *begin, const esint *end, bool omitLower);
	void clearMatrices(Builder::Request matrices, AssemblerData *data);

	void fillPermutedSparseData(double *target, const std::vector<esint> &indices, const std::vector<esint> &permutation, const std::vector<double> &values);

	std::vector<esint> _dirichletMap, _dirichletPermutation;
};

}



#endif /* SRC_PHYSICS_COMPOSER_COMPOSER_H_ */
