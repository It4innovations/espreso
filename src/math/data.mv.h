
#ifndef SRC_WRAPPERS_MATH_DATAMV_H_
#define SRC_WRAPPERS_MATH_DATAMV_H_

#include "matrix.csr.h"
#include "vector.dense.h"

#include <vector>

namespace espreso {

class MatrixCSRDistributed;
class VectorDenseDistributed;

struct DataMV {
	std::vector<std::vector<double> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > send, recv;
	std::vector<int> neighbors;
	esint minCol;

	MatrixCSR m;
	VectorDense v;

	void uniformCombination(const DataMV *first, const DataMV *second, int nfirst, int nsecond);

	void init(const MatrixCSRDistributed *matrix);
	void apply(const MatrixCSRDistributed *matrix, const VectorDenseDistributed *in, VectorDenseDistributed *out);
};

}



#endif /* SRC_WRAPPERS_MATH_DATAMV_H_ */
