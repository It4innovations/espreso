
#ifndef SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_
#define SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_

#include "esinfo/eslog.h"
#include <vector>

namespace espreso {

template <typename T> class Vector_Dense;
template <typename T> class Matrix_CSR;
template <template<typename> typename Matrix, typename T> class Vector_Distributed;
template <template<typename> typename Matrix, typename T> class Matrix_Distributed;

template <template<typename> typename Struct, typename T> struct Data_Synchronization {

	void init(Matrix_Distributed<Struct, T> &m)
	{
		eslog::error("try to initialize not implemented synchronizer.\n");
	}

	void gatherFromUpper(Matrix_Distributed<Struct, T> &m)
	{

	}

	void init(Vector_Distributed<Struct, T> &m)
	{
		eslog::error("try to initialize not implemented synchronizer.\n");
	}

	void gatherFromUpper(Vector_Distributed<Struct, T> &m)
	{

	}
};

template <typename T>
struct Data_Synchronization<Matrix_CSR, T> {
	std::vector<std::vector<T> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rOffset;
	std::vector<esint> nOffset;
	std::vector<int> neighbors;

	void init(Matrix_Distributed<Matrix_CSR, T> &m);
	void gatherFromUpper(Matrix_Distributed<Matrix_CSR, T> &m);
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
};

}

#endif /* SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_ */
