
#ifndef SRC_ANALYSIS_BUILDER_DIRECT_SYNCHRONIZATION_H_
#define SRC_ANALYSIS_BUILDER_DIRECT_SYNCHRONIZATION_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"

#include <vector>

namespace espreso {

template <typename T> class Matrix_Distributed;
template <template<typename, typename> typename Vector, typename T> class Vector_Distributed;

template <typename T>
struct Matrix_CSR_Sync {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Matrix_Distributed<T> &m);

	void gatherFromUpper(Matrix_Distributed<T> &m);
	void scatterToUpper(Matrix_Distributed<T> &m);
	void synchronize(Matrix_Distributed<T> &m)
	{
		gatherFromUpper(m);
		scatterToUpper(m);
	}
};

template <template<typename, typename> typename Vector, typename T>
struct Vector_Sync {
	virtual ~Vector_Sync() {}

	virtual void gatherFromUpper(Vector_Distributed<Vector, T> &v) =0;
	virtual void scatterToUpper(Vector_Distributed<Vector, T> &v) =0;
	void synchronize(Vector_Distributed<Vector, T> &v)
	{
		gatherFromUpper(v);
		scatterToUpper(v);
	}
};

template <typename T>
struct Vector_Dense_Sync: Vector_Sync<Vector_Dense, T> {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Vector_Distributed<Vector_Dense, T> &v);

	void gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v);
	void scatterToUpper(Vector_Distributed<Vector_Dense, T> &v);
};

template <typename T>
struct Vector_Sparse_Sync: Vector_Sync<Vector_Sparse, T> {

	void init(Vector_Distributed<Vector_Sparse, T> &v);

	void gatherFromUpper(Vector_Distributed<Vector_Sparse, T> &v);
	void scatterToUpper(Vector_Distributed<Vector_Sparse, T> &v);
};

}

#endif /* SRC_ANALYSIS_BUILDER_DIRECT_SYNCHRONIZATION_H_ */
