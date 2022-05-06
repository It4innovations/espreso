
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

		CHOLMOD(): A(new cholmod_sparse()), b(new cholmod_dense()), L(nullptr) { }
		~CHOLMOD() { delete A; delete b; }
	} cholmod;
};

namespace math {

template <>
void initSolver(Matrix_CSR<double> &A)
{
	A._solver = new Matrix_CSR_Solver();
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_start<esint>(A._solver->cholmod.common); break;
	default:
		break; // UMPAPCK
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
		break; // UMPAPCK
	}
}

template <>
void symbolicFactorization(const Matrix_CSR<double> &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		set(A._solver->cholmod.A, A);
		_analyze<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMPAPCK
	}
}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		set(A._solver->cholmod.A, A);
		_analyze<esint>(A._solver->cholmod.L, A._solver->cholmod.A, A._solver->cholmod.common);
		break;
	default:
		break; // UMPAPCK
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
		break; // UMPAPCK
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
		break; // UMPAPCK
	}
}

template <>
void solve(const Matrix_CSR<double> &A, Vector_Dense<double> &b, Vector_Dense<double> &x)
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
		break; // UMPAPCK
	}
}

template <>
void solve(const Matrix_CSR<double> &A, Matrix_Dense<double> &b, Matrix_Dense<double> &x)
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
		break; // UMPAPCK
	}
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Vector_Dense<std::complex<double> > &b, Vector_Dense<std::complex<double> > &x)
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
		break; // UMPAPCK
	}
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Matrix_Dense<std::complex<double> > &b, Matrix_Dense<std::complex<double> > &x)
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
		break; // UMPAPCK
	}
}

template <>
void freeSolver(Matrix_CSR<double> &A)
{
	switch (A.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_free<esint>(A._solver->cholmod.L, A._solver->cholmod.common);
		_finish<esint>(A._solver->cholmod.common); break;
	default:
		break; // UMPAPCK
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
		_finish<esint>(A._solver->cholmod.common); break;
	default:
		break; // UMPAPCK
	}
	delete A._solver;
}

}
}

#endif
#endif

