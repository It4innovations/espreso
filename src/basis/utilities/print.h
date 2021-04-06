
#ifndef SRC_BASIS_UTILITIES_PRINT_H_
#define SRC_BASIS_UTILITIES_PRINT_H_

#include "esinfo/mpiinfo.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "mesh/store/contactstore.h"
#include "feti/generic/SparseMatrix.h"
#include "math/matrix.dense.h"
#include "math/matrix.dense.distributed.h"
#include "math/matrix.csr.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"
#include "math/vector.sparse.h"
#include "math/domainindices.h"
#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
//	os << "<" << p.x << " " << p.y << ">";
	os << "<" << p.x << " " << p.y << " " << p.z << ">";
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const Triangle &p) {
	os << "[" << p.p[0] << p.p[1] << p.p[2] << "]";
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorDenseDistributed &v)
{
	for (esint i = v.nhalo; i < v.size; i++) {
		os << v.vals[i] << " ";
	}
	os << std::endl;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorSparse &m)
{
	for (esint i = 0; i < m.nnz; i++) {
		os << m.indices[i] << " ";
		os << std::scientific << m.vals[i] << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorDense &v)
{
	if (dynamic_cast<const VectorDenseDistributed*>(&v) != NULL) {
		os << dynamic_cast<const VectorDenseDistributed&>(v);
		return os;
	}
	if (dynamic_cast<const VectorSparse*>(&v) != NULL) {
		os << dynamic_cast<const VectorSparse&>(v);
		return os;
	}
	for (esint i = 0; i < v.size; i++) {
		os << v.vals[i] << " ";
	}
	os << std::endl;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsDense &v)
{
	for (esint n = 0; n < v.nvectors; n++) {
		os << v[n];
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixCSRDistributed &m)
{
	if (m.nnz == 0) {
		return os << m.nrows << " " << m.ncols << " 0\n";
	}
	os << m.nrows << " " << m.ncols << " " << m.nnz - m.rows[m.nhalo] + m.rows[0] << "\n";

	for (esint r = 0; r < m.nhalo; r++) {
		for (esint c = m.rows[r]; c < m.rows[r + 1]; c++) {
			os << m.halo[r] + m.rows[0] << " ";
			os << m.cols[c - m.rows[0]] << " ";
			os << std::scientific << m.vals[c - m.rows[0]] << "\n";
		}
	}

	for (esint r = m.nhalo; r < m.nrows; r++) {
		for (esint c = m.rows[r]; c < m.rows[r + 1]; c++) {
			os << r + m.rows[0] + m.distribution[info::mpi::rank] - m.nhalo << " ";
			os << m.cols[c - m.rows[0]] << " ";
			os << std::scientific << m.vals[c - m.rows[0]] << "\n";
		}
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixCSR &m)
{
	if (m.nnz == 0) {
		return os << m.nrows << " " << m.ncols << " 0\n";
	}
	if (dynamic_cast<const MatrixCSRDistributed*>(&m) != NULL) {
		os << dynamic_cast<const MatrixCSRDistributed&>(m);
		return os;
	}
	os << m.nrows << " " << m.ncols << " " << m.nnz << "\n";

	for (esint r = 0; r < m.nrows; r++) {
		for (esint c = m.rows[r]; c < m.rows[r + 1]; c++) {
			os << r + m.rows[0] << " ";
			os << m.cols[c - m.rows[0]] << " ";
			os << std::scientific << m.vals[c - m.rows[0]] << "\n";
		}
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixDense &m)
{
	for (esint i = 0; i < m.nrows; i++) {
		for (esint j = 0; j < m.ncols; j++) {
			os << m.vals[i * m.ncols + j] << " ";
		}
		os << std::endl;
	}
	os << std::endl;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixDenseDistributed &m)
{
	for (esint i = m.nhalo; i < m.nrows; i++) {
		for (esint j = 0; j < m.ncols; j++) {
			os << m.vals[i * m.ncols + j] << " ";
		}
		os << std::endl;
	}
	os << std::endl;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const SparseMatrix &m)
{
	os << m.rows << " " << m.cols << " " << m.nnz << "\n";

	os.precision(15);

	SparseMatrix s = m;
	if (s.CSR_J_col_indices.size()) {
		s.ConvertToCOO(1);
	}
	if (s.dense_values.size()) {
		s.ConvertDenseToCSR(0);
		s.ConvertToCOO(1);
	}

	for (esint i = 0; i < s.nnz; i++) {
		os << s.I_row_indices[i] << " ";
		os << s.J_col_indices[i] << " ";
		os << std::scientific << s.V_values[i] << "\n";
	}
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">";
	return os;
}

template<typename T, typename A>
std::ostream& operator<< (std::ostream& os, const std::vector<T, A> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

template <typename TData>
std::ostream& operator<< (std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<< (std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	for(auto e = data.cbegin(); e != data.cend(); ++e) {
		os << *e;
	}
	return os;
}

inline std::ostream& operator<< (std::ostream& os, const DI &di)
{
	os << di.domain << ":" << di.index;
	return os;
}

}

#endif /* SRC_BASIS_UTILITIES_PRINT_H_ */
