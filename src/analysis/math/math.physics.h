
#ifndef SRC_ANALYSIS_MATH_MATH_PHYSICS_H_
#define SRC_ANALYSIS_MATH_MATH_PHYSICS_H_

#include "selection.h"
#include "math/math.h"

namespace espreso {
namespace math {

template <typename T> void copy(Matrix_CSR<T> &x, const Matrix_CSR<T> &y, const Selection &rows, const Selection &cols)
{
	if (rows == Selection() && cols == Selection()) {
		math::copy(x, y);
	} else {
		if (y.nrows < x.nrows || y.ncols < x.ncols) {
			for (esint r = 0; r < y.nrows; r++) {
				for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
					esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
					esint tc = cols.step * (c / cols.size) + cols.offset + c % cols.size;
					x.vals[x.rows[tr] + tc - x.rows[0]] = y.vals[y.rows[r] + c - y.rows[0]];
				}
			}
		} else {
			for (esint r = 0; r < x.nrows; r++) {
				for (esint c = 0; c < x.rows[r + 1] - x.rows[r]; c++) {
					esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
					esint tc = cols.step * (c / cols.size) + cols.offset + c % cols.size;
					x.vals[x.rows[r] + c - x.rows[0]] = y.vals[y.rows[tr] + tc - y.rows[0]];
				}
			}
		}
	}
}

template <typename T> void copy(Vector_Dense<T> &x, const Vector_Dense<T> &y, const Selection &rows)
{
	if (rows == Selection()) {
		math::copy(x, y);
	} else {
		if (y.size < x.size) {
			for (esint r = 0; r < y.size; r++) {
				esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
				x.vals[tr] = y.vals[r];
			}
		} else {
			for (esint r = 0; r < x.size; r++) {
				esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
				x.vals[r] = y.vals[tr];
			}
		}
	}
}

template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Sparse<T> &y, const Selection &rows)
{
	if (rows == Selection()) {
		math::copy(x, y);
	} else {
		if (x.nnz < y.nnz) {
			for (esint i = 0; i < x.nnz / rows.size; ++i) {
				for (esint j = 0; j < rows.size; ++j) {
					x.vals[i * rows.size + j] = y.vals[i * rows.step + rows.offset + j];
				}
			}
		} else {
			for (esint i = 0; i < y.nnz / rows.size; ++i) {
				for (esint j = 0; j < rows.size; ++j) {
					x.vals[i * rows.step + rows.offset + j] = y.vals[i * rows.size + j];
				}
			}
		}
	}
}

template <typename T> void copy(Vector_Dense<T> &x, const Vector_Sparse<T> &y, const Selection &rows)
{
	for (esint i = 0; i < y.nnz; ++i) {
		x.vals[y.indices[i]] = y.vals[i];
	}
}

template <typename T> void copy(Vector_Sparse<T> &x, const Vector_Dense<T> &y, const Selection &rows)
{
	for (esint i = 0; i < x.nnz; ++i) {
		x.vals[i] = y.vals[x.indices[i]];
	}
}

template <typename T> void add(Matrix_CSR<T> &x, const T &alpha, const Matrix_CSR<T> &y, const Selection &rows, const Selection &cols)
{
	if (rows == Selection() && cols == Selection()) {
		math::add(x, alpha, y);
	} else {
		if (y.nrows < x.nrows || y.ncols < x.ncols) {
			for (esint r = 0; r < y.nrows; r++) {
				for (esint c = 0; c < y.rows[r + 1] - y.rows[r]; c++) {
					esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
					esint tc = cols.step * (c / cols.size) + cols.offset + c % cols.size;
					x.vals[x.rows[tr] + tc - x.rows[0]] += alpha * y.vals[y.rows[r] + c - y.rows[0]];
				}
			}
		} else {
			for (esint r = 0; r < x.nrows; r++) {
				for (esint c = 0; c < x.rows[r + 1] - x.rows[r]; c++) {
					esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
					esint tc = cols.step * (c / cols.size) + cols.offset + c % cols.size;
					x.vals[x.rows[r] + c - x.rows[0]] *= alpha * y.vals[y.rows[tr] + tc - y.rows[0]];
				}
			}
		}
	}
}

template <typename T> void add(Vector_Dense<T> &x, const T &alpha, const Vector_Dense<T> &y, const Selection &rows)
{
	if (rows == Selection()) {
		math::add(x, alpha, y);
	} else {
		if (y.size < x.size) {
			for (esint r = 0; r < y.size; r++) {
				esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
				x.vals[tr] += alpha * y.vals[r];
			}
		} else {
			for (esint r = 0; r < x.size; r++) {
				esint tr = rows.step * (r / rows.size) + rows.offset + r % rows.size;
				x.vals[r] += alpha * y.vals[tr];
			}
		}
	}
}

template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Sparse<T> &y, const Selection &rows)
{
	if (rows == Selection()) {
		math::add(x, alpha, y);
	} else {
		if (x.nnz < y.nnz) {
			for (esint i = 0; i < x.nnz / rows.size; ++i) {
				for (esint j = 0; j < rows.size; ++j) {
					x.vals[i * rows.size + j] += alpha * y.vals[i * rows.step + rows.offset + j];
				}
			}
		} else {
			for (esint i = 0; i < y.nnz / rows.size; ++i) {
				for (esint j = 0; j < rows.size; ++j) {
					x.vals[i * rows.step + rows.offset + j] += alpha * y.vals[i * rows.size + j];
				}
			}
		}
	}
}

template <typename T> void add(Vector_Sparse<T> &x, const T &alpha, const Vector_Dense<T> &y, const Selection &rows)
{
	for (esint i = 0; i < x.nnz; ++i) {
		x.vals[i] += alpha * y.vals[x.indices[i]];
	}
}
template <typename T> void add(Vector_Dense<T> &x, const T &alpha, const Vector_Sparse<T> &y, const Selection &rows)
{
	for (esint i = 0; i < y.nnz; ++i) {
		x.vals[y.indices[i]] += alpha * y.vals[i];
	}
}

}
}

#endif /* SRC_ANALYSIS_MATH_MATH_PHYSICS_H_ */
