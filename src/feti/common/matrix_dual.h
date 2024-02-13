
#ifndef SRC_FETI_COMMON_MATRIX_DUAL_H_
#define SRC_FETI_COMMON_MATRIX_DUAL_H_

#include "dual_map.h"
#include "feti/feti.h"
#include "math/primitives/matrix_dense.h"

#include <vector>

namespace espreso {

template <typename T>
struct Matrix_Dual: public Matrix_Dense<T> {

	void resize(int nrhs);
	void synchronize();

	int nhalo;
	using Matrix_Dense<T>::nrows;
	using Matrix_Dense<T>::ncols;
	using Matrix_Dense<T>::vals;

    static std::vector<std::vector<T> > send, recv; // buffer
};

template <typename T> std::vector<std::vector<T> > Matrix_Dual<T>::send;
template <typename T> std::vector<std::vector<T> > Matrix_Dual<T>::recv;

}

#endif /* SRC_FETI_COMMON_MATRIX_DUAL_H_ */
