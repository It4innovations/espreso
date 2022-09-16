
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef USE_SOLVER_SUITESPARSE

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"

namespace espreso {

struct Matrix_CSR_Solver {
	struct CHOLMOD {
		cholmod_sparse *A;
		cholmod_dense *b;

		cholmod_factor *L;
		cholmod_common common;

		CHOLMOD(): A(nullptr), b(nullptr), L(nullptr) { }
		~CHOLMOD() { if (A) delete A; if (b) delete b; if (L) delete L; }
	} cholmod;
};

struct Matrix_CSC_Solver: public Matrix_CSR_Solver {};

namespace math {

const char* sparseSolver()
{
	return "SUITE SPARSE";
}

template <>
void initSolver(Matrix_CSR<double> &A)
{
	A._solver = new Matrix_CSR_Solver();
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_start<esint>(A._solver->cholmod.common); break;
	default:
		break; // UMFPACK
	}
}

template <>
void initSolver(Matrix_CSR<std::complex<double> > &A)
{
	A._solver = new Matrix_CSR_Solver();
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_start<esint>(A._solver->cholmod.common); break;
	default:
		break; // UMFPACK
	}
}

template <>
void symbolicFactorization(const Matrix_CSR<double> &A, esint fixedSuffix)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		set(A._solver->cholmod.A, A);
		_analyze<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &A, esint fixedSuffix)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		set(A._solver->cholmod.A, A);
		_analyze<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void numericalFactorization(const Matrix_CSR<double> &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.A, A);
		_factorize<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void numericalFactorization(const Matrix_CSR<std::complex<double> > &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.A, A);
		_factorize<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void solve(const Matrix_CSR<double> &A, Vector_Dense<double> &b, Vector_Dense<double> &x, VectorSparsity sparsity)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.b, b);
		cholmod_dense *_x;
		_solve<esint>(_x, A._solver->cholmod.L, A._solver->cholmod.b, A._solver->cholmod.common);
		extract(_x, A._solver->cholmod.common, x);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void solve(const Matrix_CSR<double> &A, Matrix_Dense<double> &b, Matrix_Dense<double> &x, VectorSparsity sparsity)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.b, b);
		cholmod_dense *_x;
		_solve<esint>(_x, A._solver->cholmod.L, A._solver->cholmod.b, A._solver->cholmod.common);
		extract(_x, A._solver->cholmod.common, x);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Vector_Dense<std::complex<double> > &b, Vector_Dense<std::complex<double> > &x, VectorSparsity sparsity)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.b, b);
		cholmod_dense *_x;
		_solve<esint>(_x, A._solver->cholmod.L, A._solver->cholmod.b, A._solver->cholmod.common);
		extract(_x, A._solver->cholmod.common, x);
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Matrix_Dense<std::complex<double> > &b, Matrix_Dense<std::complex<double> > &x, VectorSparsity sparsity)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(A._solver->cholmod.b, b);
		cholmod_dense *_x;
		_solve<esint>(_x, A._solver->cholmod.L, A._solver->cholmod.b, A._solver->cholmod.common);
		extract(_x, A._solver->cholmod.common, x);
		break;
	default:
		break; // UMFPACK
	}
}

bool provideSC()
{
	return false;
}

template <>
void computeSC(const Matrix_CSR<double> &A, Matrix_Dense<double> &sc)
{
	eslog::error("Implement Schur complement via SuiteSparse.\n");
}

template <>
void computeSC(const Matrix_CSR<std::complex<double> > &A, Matrix_Dense<std::complex<double> > &sc)
{
	eslog::error("Implement Schur complement via SuiteSparse.\n");
}

template <>
void freeSolver(Matrix_CSR<double> &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_free<esint>(A._solver->cholmod.L, A._solver->cholmod.common);
		_finish<esint>(A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
	delete A._solver;
}

template <>
void freeSolver(Matrix_CSR<std::complex<double> > &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_free<esint>(A._solver->cholmod.L, A._solver->cholmod.common);
		_finish<esint>(A._solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
	delete A._solver;
}

template<typename T>
static void _info(SolverInfo &info, const Matrix_CSR<T> &A)
{
	info.rows = A.nrows;
	info.nnzA = A.nnz;
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		info.nnzL = A._solver->cholmod.common.lnz;
		info.memoryL = A._solver->cholmod.common.memory_inuse;
		break;
	default:
		info.nnzL = 0;
		break; // UMFPACK
	}
}

template <>
SolverInfo getSolverInfo(const Matrix_CSR<double> &A)
{
	SolverInfo info;
	_info(info, A);
	return info;
}

template <>
SolverInfo getSolverInfo(const Matrix_CSR<std::complex<double> > &A)
{
	SolverInfo info;
	_info(info, A);
	return info;
}

bool provideFactors()
{
	return true;
}

template <>
void getFactors(const Matrix_CSR<double> &A, Matrix_CSC<double> &L, Matrix_CSC<double> &U, Vector_Dense<int> &p)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		L._solver = new Matrix_CSC_Solver();
		L._solver->cholmod.common = A._solver->cholmod.common;
		cholmod_factor *copy;
		_copyFactor<esint>(A._solver->cholmod.L, copy, A._solver->cholmod.common);
		_factorToSparse<esint>(copy, L._solver->cholmod.A, L._solver->cholmod.common);
		_free<esint>(copy, A._solver->cholmod.common);
		L.nrows = L._solver->cholmod.A->nrow;
		L.ncols = L._solver->cholmod.A->ncol;
		L.nnz = L._solver->cholmod.A->nzmax;
		L.rows = (esint*)L._solver->cholmod.A->p;
		L.cols = (esint*)L._solver->cholmod.A->i;
		L.vals = (double*)L._solver->cholmod.A->x;
		L.type = A.type;
		L.shape = Matrix_Shape::LOWER;
//		p.size = L.nrows;
//		p.vals = (esint*)A._solver->cholmod.L->Perm;
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void getFactors(const Matrix_CSR<std::complex<double> > &A, Matrix_CSC<std::complex<double> > &L, Matrix_CSC<std::complex<double> > &U, Vector_Dense<int> &p)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		L._solver = new Matrix_CSC_Solver();
//		L._solver->cholmod.common = A._solver->cholmod.common;

		_factorToSparse<esint>(A._solver->cholmod.L, L._solver->cholmod.A, A._solver->cholmod.common);
		L.nrows = L._solver->cholmod.A->nrow;
		L.ncols = L._solver->cholmod.A->ncol;
		L.nnz = L._solver->cholmod.A->nzmax;
		L.rows = (esint*)L._solver->cholmod.A->p;
		L.cols = (esint*)L._solver->cholmod.A->i;
		L.vals = (std::complex<double>*)L._solver->cholmod.A->z;
		L.type = A.type;
		L.shape = Matrix_Shape::LOWER;
		p.size = L.nrows;
		p.vals = (esint*)A._solver->cholmod.L->Perm;
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void freeFactor(Matrix_CSC<double> &L)
{
	switch (L.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_free<esint>(L._solver->cholmod.A, L._solver->cholmod.common);
		delete L._solver;
		break;
	default:
		break; // UMFPACK
	}
}

template <>
void freeFactor(Matrix_CSC<std::complex<double> > &L)
{
	switch (L.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_free<esint>(L._solver->cholmod.A, L._solver->cholmod.common);
		delete L._solver;
		break;
	default:
		break; // UMFPACK
	}
}

}
}

#endif
#endif

