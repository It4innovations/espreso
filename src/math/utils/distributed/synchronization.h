
#ifndef SRC_MATH2_UTILS_DISTRIBUTED_SYNCHRONIZATION_H_
#define SRC_MATH2_UTILS_DISTRIBUTED_SYNCHRONIZATION_H_

#include "esinfo/eslog.h"

#include "math/primitives/matrix_csr.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "wrappers/mpi/communication.h"
#include "mesh/store/nodestore.h"

#include <vector>

namespace espreso {

template <template<typename, typename, template<typename> typename> typename Vector, typename T> class Vector_Distributed;
template <template<typename, typename, template<typename> typename> typename Matrix, typename T> class Matrix_Distributed;

template <template<typename, typename, template<typename> typename> typename Struct, typename T> struct Data_Synchronization { };

template <typename T>
struct Data_Synchronization<Matrix_CSR, T> {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Matrix_Distributed<Matrix_CSR, T> &m);
	void gatherFromUpper(Matrix_Distributed<Matrix_CSR, T> &m);
	void scatterToUpper(Matrix_Distributed<Matrix_CSR, T> &m);

	void synchronize(Matrix_Distributed<Matrix_CSR, T> &m)
	{
		gatherFromUpper(m);
		scatterToUpper(m);
	}
};

template <typename T>
struct Data_Synchronization<Vector_Dense, T> {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Vector_Distributed<Vector_Dense, T> &v);
	void gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v);
	void scatterToUpper(Vector_Distributed<Vector_Dense, T> &v);

	void synchronize(Vector_Distributed<Vector_Dense, T> &v)
	{
		gatherFromUpper(v);
		scatterToUpper(v);
	}
};

template <typename T>
struct Data_Synchronization<Vector_Sparse, T> {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Vector_Distributed<Vector_Sparse, T> &v) { eslog::error("try to initialize not implemented synchronizer.\n"); }
	void gatherFromUpper(Vector_Distributed<Vector_Sparse, T> &v) {}
	void scatterToUpper(Vector_Distributed<Vector_Sparse, T> &v) {}

	void synchronize(Vector_Distributed<Vector_Sparse, T> &v)
	{
		gatherFromUpper(v);
		scatterToUpper(v);
	}
};

}

#endif /* SRC_MATH2_UTILS_DISTRIBUTED_SYNCHRONIZATION_H_ */
