
#include "math/math.h"
#include "math/primitives/vector_dense.h"
#include "math/feti/lmap.h"
#include "math/feti/vector_dual.h"
#include "math/feti/vector_kernel.h"
#include "math/generalization/vector_distributed.h"
#include "math/generalization/vector_feti.h"
#include "math/generalization/matrix_distributed.h"
#include "math/generalization/matrix_feti.h"
#include "esinfo/meshinfo.h"

#include <fstream>
#include <iomanip>

namespace espreso {
namespace math {

template <typename T>
void _store(const Vector_Dense<T> &x, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os.precision(15);
	os << std::showpos;

	for (esint i = 0; i < x.size; i++) {
		os << std::setw(25) << std::scientific << x.vals[i] << "\n";
	}
}

template <> void store(const Vector_Dense<double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_Dense<std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Vector_Sparse<T> &x, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os.precision(15);
	os << std::showpos;

	for (esint i = 0; i < x.nnz; i++) {
		os << std::setw(6) << x.indices[i] << " ";
		os << std::setw(25) << std::scientific << x.vals[i] << "\n";
	}
}

template <> void store(const Vector_Sparse<double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_Sparse<std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Vector_Distributed<Vector_Dense, T> &x, const char* file)
{
	Vector_Dense<T> _x;
	_x.size = x.cluster.size - x.distribution->halo.size();
	_x.vals = x.cluster.vals + x.distribution->halo.size();
	store(_x, file);
}

template <> void store(const Vector_Distributed<Vector_Dense, double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_Distributed<Vector_Dense, std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Vector_Distributed<Vector_Sparse, T> &x, const char* file)
{
	Vector_Sparse<T> _x;
	_x.size = x.cluster.size - x.distribution->halo.size();
	_x.nnz = x.cluster.nnz - x.distribution->halo.size();
	_x.indices = x.cluster.indices + x.distribution->halo.size();
	_x.vals = x.cluster.vals + x.distribution->halo.size();
	store(_x, file);
}

template <> void store(const Vector_Distributed<Vector_Sparse, double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_Distributed<Vector_Sparse, std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Vector_FETI<Vector_Dense, T> &x, const char* file)
{
	for (size_t d = 0; d < x.domains.size(); ++d) {
		store(x.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <> void store(const Vector_FETI<Vector_Dense, double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_FETI<Vector_Dense, std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Vector_FETI<Vector_Sparse, T> &x, const char* file)
{
	for (size_t d = 0; d < x.domains.size(); ++d) {
		store(x.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <> void store(const Vector_FETI<Vector_Sparse, double> &x, const char* file) { _store(x, file); }
template <> void store(const Vector_FETI<Vector_Sparse, std::complex<double> > &x, const char* file) { _store(x, file); }

template <typename T>
void _store(const Matrix_Dense<T> &A, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os << std::setw(6) << A.nrows << " " << std::setw(6) << A.ncols << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint r = 0; r < A.nrows; r++) {
		for (esint c = 0; c < A.ncols; c++) {
			os << std::setw(25) << std::scientific << A.vals[r * A.ncols + c];
		}
		os << "\n";
	}
}

template <> void store(const Matrix_Dense<double> &A, const char* file) { _store(A, file); }
template <> void store(const Matrix_Dense<std::complex<double> > &A, const char* file) { _store(A, file); }

template <typename T>
void _store(const Matrix_CSR<T> &A, const char* file, esint offset = 0)
{
	std::ofstream os(std::string(file) + ".txt");
	os << std::setw(6) << A.nrows << " " << std::setw(6) << A.ncols << " " << std::setw(6) << A.nnz << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint r = 0; r < A.nrows; r++) {
		for (esint c = A.rows[r]; c < A.rows[r + 1]; c++) {
			os << std::setw(6) << offset + r + A.rows[0] << " ";
			os << std::setw(6) << A.cols[c - A.rows[0]] << " ";
			os << std::setw(25) << std::scientific << A.vals[c - A.rows[0]] << "\n";
		}
	}
}

template <> void store(const Matrix_CSR<double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_CSR<std::complex<double>> &A, const char* file) { _store<std::complex<double>>(A, file); }

template <typename T>
void _store(const Matrix_IJV<T> &A, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os << std::setw(6) << A.nrows << " " << std::setw(6) << A.ncols << " " << std::setw(6) << A.nnz << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint i = 0; i < A.nnz; i++) {
		os << std::setw(6) << A.rows[i]<< " ";
		os << std::setw(6) << A.cols[i] << " ";
		os << std::setw(25) << std::scientific << A.vals[i] << "\n";
	}
}

template <> void store(const Matrix_IJV<double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_IJV<std::complex<double>> &A, const char* file) { _store<std::complex<double>>(A, file); }

template <typename T>
void _store(const Matrix_Distributed<Matrix_CSR, T> &A, const char* file)
{
	Matrix_CSR<T> _A;
	_A.nrows = A.cluster.nrows - A.distribution->halo.size();
	_A.ncols = A.cluster.ncols;
	_A.nnz = A.cluster.nnz - (A.cluster.rows[A.distribution->halo.size()] - Indexing::CSR);
	_A.rows = A.cluster.rows + A.distribution->halo.size();
	_A.cols = A.cluster.cols;
	_A.vals = A.cluster.vals;
	_store(_A, file, A.distribution->begin);
}

template <> void store(const Matrix_Distributed<Matrix_CSR, double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_Distributed<Matrix_CSR, std::complex<double> > &A, const char* file) { _store<std::complex<double>>(A, file); }

template <typename T>
void _store(const Matrix_FETI<Matrix_Dense, T> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <> void store(const Matrix_FETI<Matrix_Dense, double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_FETI<Matrix_Dense, std::complex<double> > &A, const char* file) { _store<std::complex<double>>(A, file); }

template <typename T>
void _store(const Matrix_FETI<Matrix_CSR, T> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <> void store(const Matrix_FETI<Matrix_CSR, double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_FETI<Matrix_CSR, std::complex<double> > &A, const char* file) { _store<std::complex<double>>(A, file); }

template <typename T>
void _store(const Matrix_FETI<Matrix_IJV, T> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <> void store(const Matrix_FETI<Matrix_IJV, double> &A, const char* file) { _store<double>(A, file); }
template <> void store(const Matrix_FETI<Matrix_IJV, std::complex<double> > &A, const char* file) { _store<std::complex<double>>(A, file); }

template <>
void store(const std::vector<esint> &v, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");

	for (size_t i = 0; i < v.size(); i++) {
		os << std::setw(6) << v[i] << "\n";
	}
}

template <>
void store(const std::vector<std::vector<esint> > &v, const char* file)
{
	for (size_t d = 0; d < v.size(); ++d) {
		store(v[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <>
void store(const std::vector<LMAP> &lmap, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	for (auto l = lmap.cbegin(); l != lmap.cend(); ++l) {
		if (l->neigh < 0) {
			os << std::setw(4) << l->offset << ": " << l->from << "->" << l->to << "\n";
		} else {
			os << std::setw(4) << l->offset << ": " << l->from << "->" << l->to << "::" << info::mesh->neighbors[l->neigh] << "\n";
		}
	}
}

template <>
void store(const Vector_Dual<double> &x, const char* file)
{
	math::store(static_cast<espreso::Vector_Dense<double> >(x), file);
}

template <>
void store(const Vector_Dual<std::complex<double> > &x, const char* file)
{
	math::store(static_cast<espreso::Vector_Dense<std::complex<double> > >(x), file);
}

template <>
void store(const Vector_Kernel<double> &x, const char* file)
{
	math::store(static_cast<espreso::Vector_Dense<double> >(x), file);
}

template <>
void store(const Vector_Kernel<std::complex<double> > &x, const char* file)
{
	math::store(static_cast<espreso::Vector_Dense<std::complex<double> > >(x), file);
}

}
}
