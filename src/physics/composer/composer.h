
#ifndef SRC_PHYSICS_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_COMPOSER_H_

#include "physics/system/builder/builder.h"
#include <cstddef>
#include <vector>

namespace espreso {

class Kernel;
struct AssemblerData;
struct IJ;
class VectorSparse;

class Composer {

public:
	Composer(Kernel *kernel);

	virtual void init() = 0;
	virtual void assemble(const Builder &builder) = 0;

	virtual ~Composer();

	Kernel *kernel;

protected:
	esint getMatrixSize(esint size, bool omitLower);
	void insertKPattern(IJ *target, const esint *begin, const esint *end, bool omitLower);
	void clearMatrices(Builder::Request matrices, AssemblerData *data);

	void fillPermutedSparseData(double *target, const std::vector<esint> &indices, const std::vector<esint> &permutation, const std::vector<double> &values);

	std::vector<esint> _dirichletMap, _dirichletPermutation;
};

}



#endif /* SRC_PHYSICS_COMPOSER_COMPOSER_H_ */
