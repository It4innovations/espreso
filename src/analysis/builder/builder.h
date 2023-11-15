
#ifndef SRC_ANALYSIS_BUILDER_BUILDER_H_
#define SRC_ANALYSIS_BUILDER_BUILDER_H_

#include "math/physics/matrix_base.h"
#include "math/physics/vector_base.h"

#include <map>
#include <vector>

namespace espreso {

template <typename T>
struct SparseMatrixBuilder {

	SparseMatrixBuilder(int dofs): dofs(dofs) {}
	virtual ~SparseMatrixBuilder() {}

	virtual void fillDirichlet(Vector_Base<T> *v) =0;
	virtual void fillVector(Vector_Base<T> *v) =0;
	virtual void fillMatrix(Matrix_Base<T> *m, Matrix_Type type, Matrix_Shape shape) =0;

	virtual void fillDirichletMap(Vector_Base<T> *v) =0;
	virtual void fillVectorMap(Vector_Base<T> *v) =0;
	virtual void fillMatrixMap(Matrix_Base<T> *m) =0;

	int dofs;
};

}

#endif /* SRC_ANALYSIS_BUILDER_BUILDER_H_ */