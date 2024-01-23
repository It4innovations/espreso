
#ifndef SRC_MATH_FETI_MATRIX_DUAL_ORTHOGONAL_H_
#define SRC_MATH_FETI_MATRIX_DUAL_ORTHOGONAL_H_

#include "vector_dual.h"
#include "esinfo/envinfo.h"
#include "math/math.h"
#include "math/primitives/matrix_dense.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <memory>

namespace espreso {

template <typename T>
struct Matrix_Dual_Orthogonal: public Matrix_Dense<T, int> {

	enum: size_t {
		initial_space = 32
	};

	void resize()
	{
		// align matrix values ??
		Matrix_Dense<T>::resize(initial_space, Vector_Dual<T>::localSize);
		Matrix_Dense<T>::nrows = 0;
		Matrix_Dense<T>::slice({}, { Vector_Dual<T>::nhalo, Vector_Dual<T>::localSize });
	}

	void next(Vector_Dual<T> &v)
	{
		if (Matrix_Dense<T>::nrows == Matrix_Dense<T>::_allocated.nrows) {
			Matrix_Dense<T> _m;
			_m.resize(2 * Matrix_Dense<T>::nrows, Matrix_Dense<T>::ncols);
			_m.nrows = Matrix_Dense<T>::nrows;
			Matrix_Dense<T>::swap(*this, _m);
			Matrix_Dense<T>::swap(Matrix_Dense<T>::_allocated, _m.allocated());
			memcpy(Matrix_Dense<T>::vals, _m.vals, sizeof(T) * _m.nrows * _m.ncols);
		}
		v.vals = Matrix_Dense<T>::vals + Matrix_Dense<T>::ncols * Matrix_Dense<T>::nrows++;
		v.size = Vector_Dual<T>::localSize;
	}

	void apply(const Vector_Dual<T> &x, Vector_Dense<T> &y)
	{
		Vector_Dense<T> _x;
		_x.size = x.size - Vector_Dual<T>::nhalo;
		_x.vals = x.vals + Vector_Dual<T>::nhalo;
		math::blas::apply(y, T{1}, static_cast<Matrix_Dense<T>&>(*this), T{0}, _x);
		Communication::allReduce(y.vals, nullptr, y.size, MPITools::getType<T>().mpitype, MPI_SUM);
	}

	void applyT(const Vector_Dense<T> &x, Vector_Dual<T> &y)
	{
		esint size = x.size;
		Slice slice(0, Vector_Dual<T>::localSize);
		std::swap(Matrix_Dense<T>::nrows, size);
		std::swap(Matrix_Dense<T>::submatrix[1], slice);
		math::blas::applyT(y, T{1}, static_cast<Matrix_Dense<T>&>(*this), T{0}, x);
		std::swap(Matrix_Dense<T>::nrows, size);
		std::swap(Matrix_Dense<T>::submatrix[1], slice);
	}
};

}

#endif /* SRC_MATH_FETI_MATRIX_DUAL_ORTHOGONAL_H_ */
