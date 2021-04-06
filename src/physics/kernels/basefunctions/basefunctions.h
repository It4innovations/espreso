
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_BASEFUNCTIONS_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_BASEFUNCTIONS_H_

#include "mesh/element.h"

#include <cstddef>
#include <vector>

namespace espreso {

class MatrixDense;

struct BaseFunctions: public Element {

	static void setBaseFunctions();
	static void created(Element &e);

	virtual ~BaseFunctions();

	virtual void computeReferenceCoords(const MatrixDense & vertices, const MatrixDense & points, MatrixDense & result) {}
	void recomputeDetJ( MatrixDense& coords, MatrixDense& resdetJ, MatrixDense* points = NULL);
	void recomputeDetJN(MatrixDense& coords, MatrixDense& resdetJ, MatrixDense& resN, MatrixDense& refPoints);


	virtual void setGaussPointsForOrder(int order) =0;

};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_BASEFUNCTIONS_H_ */
