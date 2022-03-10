
#ifndef SRC_AXFETI_PROJECTOR_ORTHOGONAL_H_
#define SRC_AXFETI_PROJECTOR_ORTHOGONAL_H_

#include "projector.h"

#include <unordered_map>

namespace espreso {

// G = R'B'

// 1. per process R'B' -> rows of G
// 2. Gather G
// 3. compute GG'
// 4a. Factorization of GG'
// 4b. Explicit inv(GG')

template <typename T>
class Orthogonal: public Projector<T> {
public:
	Orthogonal(AX_FETI<T> *feti);

	void info();
	void update();

	Matrix_CSR<T> nnG, localG, GGt;
	Matrix_Dense<T> invGGt;
	esint nnGPreRows;
	std::unordered_map<esint, esint> nnGMap; // global column to nnG row
	std::vector<std::vector<esint> > sIndices, rIndices;
	std::vector<std::vector<T> > sBuffer, rBuffer;
};

}


#endif /* SRC_AXFETI_PROJECTOR_ORTHOGONAL_H_ */
