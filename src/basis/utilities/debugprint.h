
#ifndef SRC_BASIS_UTILITIES_DEBUGPRINT_H_
#define SRC_BASIS_UTILITIES_DEBUGPRINT_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "feti/generic/SparseMatrix.h"
#include "math/matrix.dense.h"
#include "math/matrix.dense.distributed.h"
#include "math/matrix.csr.h"
#include "math/matrix.csr.distributed.h"
#include "math/matrix.ijv.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"
#include "math/vector.dense.feti.h"
#include "math/vector.sparse.h"
#include "math/vector.sparse.feti.h"
#include <ostream>
#include <vector>
#include <iomanip>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
	os << "<" << p.x << " " << p.y << " " << p.z << ">\n";
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const SparseMatrix &m)
{
	os << std::setw(6) << m.rows << " " << m.cols << " " << m.nnz << "\n";

	SparseMatrix s = m;
	if (s.CSR_J_col_indices.size()) {
		s.ConvertToCOO(1);
	}
	if (s.dense_values.size()) {
		s.ConvertDenseToCSR(0);
		s.ConvertToCOO(1);
	}

	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < s.nnz; i++) {
		os << std::setw(6) << s.I_row_indices[i] << " ";
		os << std::setw(6) << s.J_col_indices[i] << " ";
		os << std::setw(25) << std::scientific << s.V_values[i] << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorDenseDistributed &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = v.nhalo; i < v.size; i++) {
		os << std::setw(25) << std::scientific << v.vals[i] << std::endl;
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsDenseDistributed &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = v[0].nhalo; i < v[0].size; i++) {
		os << std::setw(25) << std::scientific << v[0].vals[i];
		for (esint n = 1; n < v.nvectors; n++) {
			os << " " << std::setw(25) << std::scientific << v[n].vals[i];
		}
		os << std::endl;
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsDenseFETI::Domain &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v.data->at(0)->at(v.domain)->size; i++) {
		os << std::setw(25) << std::scientific << v.data->at(0)->at(v.domain)->vals[i];
		for (esint n = 1; n < v.data->nvectors; n++) {
			os << " " << std::setw(25) << std::scientific << v.data->at(n)->at(v.domain)->vals[i];
		}
		os << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorSparse &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v.nnz; i++) {
		os << std::setw(6) << v.indices[i] << " ";
		os << std::setw(25) << std::scientific << v.vals[i] << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsSparse &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v[0].nnz; i++) {
		os << std::setw(6) << v[0].indices[i] << " ";
		os << std::setw(25) << std::scientific << v[0].vals[i];
		for (esint n = 1; n < v.nvectors; n++) {
			os << " " << std::setw(25) << std::scientific << v[n].vals[i];
		}
		os << std::endl;
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsSparseFETI::Domain &v)
{
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v.data->at(0)->at(v.domain)->nnz; i++) {
		os << std::setw(6) << v.data->at(0)->at(v.domain)->indices[i] << " ";
		os << std::setw(25) << std::scientific << v.data->at(0)->at(v.domain)->vals[i];
		for (esint n = 1; n < v.data->nvectors; n++) {
			os << " " << std::setw(25) << std::scientific << v.data->at(n)->at(v.domain)->vals[i];
		}
		os << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorDense &v)
{
	if (dynamic_cast<const VectorDenseDistributed*>(&v) != NULL) {
		os << dynamic_cast<const VectorDenseDistributed&>(v);
		return os;
	}
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v.size; i++) {
		os << std::setw(25) << std::scientific << v.vals[i] << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const VectorsDense &v)
{
	if (dynamic_cast<const VectorsDenseDistributed*>(&v) != NULL) {
		os << dynamic_cast<const VectorsDenseDistributed&>(v);
		return os;
	}
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < v[0].size; i++) {
		os << std::setw(25) << std::scientific << v[0].vals[i] << "\n";
		for (esint n = 1; n < v.nvectors; n++) {
			os << " " << std::setw(25) << std::scientific << v[n].vals[i];
		}
		os << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixCSRDistributed &m)
{
	if (m.nnz == 0) {
		return os << m.nrows << " " << m.ncols << " 0\n";
	}
	os << std::setw(6) << m.nrows << " " << m.ncols << " " << m.nnz - m.rows[m.nhalo] + m.rows[0] << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint r = m.nhalo; r < m.nrows; r++) {
		for (esint c = m.rows[r]; c < m.rows[r + 1]; c++) {
			os << std::setw(6) << r + m.rows[0] << " ";
			os << std::setw(6) << m.cols[c - m.rows[0]] << " ";
			os << std::setw(25) << std::scientific << m.vals[c - m.rows[0]] << "\n";
		}
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixCSR &m)
{
	if (dynamic_cast<const MatrixCSRDistributed*>(&m) != NULL) {
		os << dynamic_cast<const MatrixCSRDistributed&>(m);
		return os;
	}
	os << std::setw(6) << m.nrows << " " << m.ncols << " " << m.nnz << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint r = 0; r < m.nrows; r++) {
		for (esint c = m.rows[r]; c < m.rows[r + 1]; c++) {
			os << std::setw(6) << r + m.rows[0] << " ";
			os << std::setw(6) << m.cols[c - m.rows[0]] << " ";
			os << std::setw(25) << std::scientific << m.vals[c - m.rows[0]] << "\n";
		}
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixIJV &m)
{
	if (dynamic_cast<const MatrixCSRDistributed*>(&m) != NULL) {
		os << dynamic_cast<const MatrixCSRDistributed&>(m);
		return os;
	}
	os << std::setw(6) << m.nrows << " " << m.ncols << " " << m.nnz << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < m.nnz; i++) {
		os << std::setw(6) << m.rows[i] << " ";
		os << std::setw(6) << m.cols[i] << " ";
		os << std::setw(25) << std::scientific << m.vals[i] << "\n";
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixDense &m)
{
	if (dynamic_cast<const MatrixCSR*>(&m) != NULL) {
		os << dynamic_cast<const MatrixCSR&>(m);
		return os;
	}
	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < m.nrows; i++) {
		for (esint j = 0; j < m.ncols; j++) {
			os << std::setw(25) << std::scientific << m.vals[i * m.ncols + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const MatrixDenseDistributed &m)
{
	if (dynamic_cast<const MatrixCSR*>(&m) != NULL) {
		os << dynamic_cast<const MatrixCSR&>(m);
		return os;
	}
	os.precision(15);
	os << std::showpos;
	for (esint i = m.nhalo; i < m.nrows; i++) {
		for (esint j = 0; j < m.ncols; j++) {
			os << std::setw(25) << std::scientific << m.vals[i * m.ncols + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">\n";
	return os;
}

template<typename T, typename TAlloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, TAlloc> &v)
{
	os.precision(15);
	os << std::showpos;
	for(size_t i = 0; i < v.size(); ++i) {
		os << std::setw(25) << std::scientific << v[i] << "\n";
	}
	return os;
}

template <typename TData>
std::ostream& operator<<(std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]\n";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<<(std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	size_t i = 0;
	for(auto e = data.cbegin(); e != data.cend(); ++e, ++i) {
		os << i << ": " << *e;
	}
	return os;
}

inline std::ostream& operator<< (std::ostream& os, const DI &di)
{
	os << di.domain << ":" << di.index;
	return os;
}

}



#endif /* SRC_BASIS_UTILITIES_DEBUGPRINT_H_ */
