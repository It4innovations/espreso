
#ifndef SRC_WRAPPERS_MATH_DATASYNCHRONIZATION_H_
#define SRC_WRAPPERS_MATH_DATASYNCHRONIZATION_H_

#include <vector>

namespace espreso {

class MatrixCSRDistributed;
class VectorDenseDistributed;
class VectorSparseDistributed;

struct DataSynchronization {
	std::vector<std::vector<double> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > sIndices, rIndices;
	std::vector<int> neighbors;

	void uniformCombination(const DataSynchronization *first, const DataSynchronization *second, int nfirst, int nsecond);

	void init(const MatrixCSRDistributed *matrix);
	void init(const VectorDenseDistributed *vector);
	void init(const VectorSparseDistributed *vector);

	void gatherFromUpper(const MatrixCSRDistributed *matrix);
	void gatherFromUpper(const VectorDenseDistributed *vector);
	void gatherFromUpper(const VectorSparseDistributed *vector);

	void scatterToUpper(const MatrixCSRDistributed *matrix);
	void scatterToUpper(const VectorDenseDistributed *vector);
	void scatterToUpper(const VectorSparseDistributed *vector);
};

}


#endif /* SRC_WRAPPERS_MATH_DATASYNCHRONIZATION_H_ */
