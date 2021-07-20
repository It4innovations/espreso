
#include "math2/math2.h"
#include "math2/primitives/vector_dense.h"
#include "math2/generalization/vector_distributed.h"
#include "math2/generalization/vector_feti.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/generalization/matrix_feti.h"

#include <fstream>
#include <iomanip>

namespace espreso {
namespace math {

template <>
void store(Vector_Dense<double> &x, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os.precision(15);
	os << std::showpos;

	for (esint i = 0; i < x.size; i++) {
		os << std::setw(25) << std::scientific << x.vals[i] << "\n";
	}
}

template <>
void store(Vector_Sparse<double> &x, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os.precision(15);
	os << std::showpos;

	for (esint i = 0; i < x.nnz; i++) {
		os << std::setw(6) << x.indices[i] << " ";
		os << std::setw(25) << std::scientific << x.vals[i] << "\n";
	}
}

template <>
void store(Vector_Distributed<Vector_Dense, double> &x, const char* file)
{
	Vector_Dense<double> _x;
	_x.size = x.cluster.size - x.distribution.halo.size();
	_x.vals = x.cluster.vals + x.distribution.halo.size();
	store(_x, file);
}

template <>
void store(Vector_FETI<Vector_Dense, double> &x, const char* file)
{
	for (size_t d = 0; d < x.domains.size(); ++d) {
		store(x.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <>
void store(Matrix_Dense<double> &A, const char* file)
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

template <>
void store(Matrix_CSR<double> &A, const char* file)
{
	std::ofstream os(std::string(file) + ".txt");
	os << std::setw(6) << A.nrows << " " << std::setw(6) << A.ncols << " " << std::setw(6) << A.nnz << "\n";

	os.precision(15);
	os << std::showpos;
	for (esint r = 0; r < A.nrows; r++) {
		for (esint c = A.rows[r]; c < A.rows[r + 1]; c++) {
			os << std::setw(6) << r + _Matrix_CSR_Pattern::Indexing << " ";
			os << std::setw(6) << A.cols[c - _Matrix_CSR_Pattern::Indexing] << " ";
			os << std::setw(25) << std::scientific << A.vals[c - _Matrix_CSR_Pattern::Indexing] << "\n";
		}
	}
}

template <>
void store(Matrix_IJV<double> &A, const char* file)
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

template <>
void store(Matrix_Distributed<Matrix_CSR, double> &A, const char* file)
{
	Matrix_CSR<double> _A;
	esint prefix = A.cluster.rows[A.distribution.halo.size()] - _Matrix_CSR_Pattern::Indexing;
	_A.nrows = A.cluster.nrows - A.distribution.halo.size();
	_A.ncols = A.cluster.ncols;
	_A.nnz = A.cluster.nnz - prefix;
	_A.rows = A.cluster.rows + A.distribution.halo.size();
	_A.cols = A.cluster.cols + prefix;
	_A.vals = A.cluster.vals + prefix;
	store(_A, file);
}

template <>
void store(Matrix_FETI<Matrix_Dense, double> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <>
void store(Matrix_FETI<Matrix_CSR, double> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

template <>
void store(Matrix_FETI<Matrix_IJV, double> &A, const char* file)
{
	for (size_t d = 0; d < A.domains.size(); ++d) {
		store(A.domains[d], (std::string(file) + std::to_string(d)).c_str());
	}
}

}
}
