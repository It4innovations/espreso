
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_SUITESPARSE
#ifdef USE_SOLVER_SUITESPARSE

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"
#include "math/wrappers/math.spblas.h"

namespace espreso {

template class DirectSolver<double, Matrix_CSR>;
template class DirectSolver<std::complex<double>, Matrix_CSR>;

struct Matrix_Solver_External_Representation {
	struct CHOLMOD {
		cholmod_sparse *A;
		cholmod_dense *b;

		cholmod_factor *L;
		cholmod_common common;

		CHOLMOD(): A(nullptr), b(nullptr), L(nullptr) { }
		~CHOLMOD() { if (A) delete A; if (b) delete b; if (L) delete L; }
	} cholmod;
};

template <typename T, template <typename> class Matrix>
const char* DirectSolver<T, Matrix>::name()
{
	return "SUITE SPARSE";
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideFactors()
{
	return true;
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideSC()
{
	// manually computed
	return true;
}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver()
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{

}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::~DirectSolver()
{
	if (_solver) {
		switch (matrix->type) {
		case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
			_free<esint>(_solver->cholmod.L, _solver->cholmod.common);
			_finish<esint>(_solver->cholmod.common);
			break;
		default:
			break; // UMFPACK
		}
		delete _solver;
	}
}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver(const Matrix<T> &a)
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	commit(a);
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::commit(const Matrix<T> &a)
{
	matrix = &a;
	if (_solver) { delete _solver; }
	_solver = new Matrix_Solver_External_Representation();
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_start<esint>(_solver->cholmod.common); break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::commit(SpBLAS<T, Matrix> &spblas)
{
	if (_solver) { delete _solver; }
	_solver = new Matrix_Solver_External_Representation();
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		_start<esint>(_solver->cholmod.common); break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::symbolicFactorization(int fixedSuffix)
{
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		setSymmetric(_solver->cholmod.A, *matrix);
		_analyze<esint>(_solver->cholmod.L, _solver->cholmod.A, _solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
	rows = matrix->nrows;
	nnzA = matrix->nnz;
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		nnzL = _solver->cholmod.common.lnz;
		memoryL = _solver->cholmod.common.memory_inuse;
		break;
	default:
		nnzL = 0;
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::numericalFactorization()
{
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		updateSymmetric(_solver->cholmod.A, *matrix);
		_factorize<esint>(_solver->cholmod.L, _solver->cholmod.A, _solver->cholmod.common);
		break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Vector_Dense<T> &rhs, Vector_Dense<T> &solution, int sparsity)
{
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(_solver->cholmod.b, rhs);
		cholmod_dense *_x;
		_solve<esint>(_x, _solver->cholmod.L, _solver->cholmod.b, _solver->cholmod.common);
		extract(_x, _solver->cholmod.common, solution);
		break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution, int sparsity)
{
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
		update(_solver->cholmod.b, rhs);
		cholmod_dense *_x;
		_solve<esint>(_x, _solver->cholmod.L, _solver->cholmod.b, _solver->cholmod.common);
		extract(_x, _solver->cholmod.common, solution);
		break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::getFactors(Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p)
{
	switch (matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE: {
		auto Lsolver = Matrix_Solver_External_Representation();
		Lsolver.cholmod.common = _solver->cholmod.common;
		cholmod_factor *copy;
		_copyFactor<esint>(_solver->cholmod.L, copy, _solver->cholmod.common);
		_factorToSparse<esint>(copy, Lsolver.cholmod.A, Lsolver.cholmod.common);
		_free<esint>(copy, _solver->cholmod.common);
		L.nrows = Lsolver.cholmod.A->nrow;
		L.ncols = Lsolver.cholmod.A->ncol;
		L.nnz = Lsolver.cholmod.A->nzmax;
		L.rows = (esint*)Lsolver.cholmod.A->p;
		L.cols = (esint*)Lsolver.cholmod.A->i;
		L.vals = (T*)Lsolver.cholmod.A->x;
		L.type = matrix->type;
		L.shape = Matrix_Shape::LOWER;
//		p.size = L.nrows;
//		p.vals = (esint*)A._solver->cholmod.L->Perm;
	} break;
	default:
		break; // UMFPACK
	}
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::getSC(Matrix_Dense<T> &sc)
{
	// computes the schur complement S = A22 - A21 * A11^{-1} * A12, where A = [A11, A12; A21, A22]
	switch(matrix->type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
	{
		if(matrix->shape != Matrix_Shape::UPPER) {
			eslog::error("Implement Schur complement for non-upper csr matrices.\n");
		}

		esint size_sc = sc.nrows;
		esint size = matrix->nrows;
		esint size_A11 = size - size_sc;

		Matrix_CSR<T> A11_sp;
		Matrix_CSR<T> A21t_sp; // = A12c_sp
		Matrix_Dense<T> A22t_dn;
		Matrix_Dense<T> A12t_dn;
		SpBLAS<T, Matrix_CSR>::submatrix(*matrix, A11_sp, 0, size_A11, 0, size_A11);
		SpBLAS<T, Matrix_CSR>::submatrix(*matrix, A21t_sp, 0, size_A11, size_A11, size, false, true); // = A12c_sp
		SpBLAS<T, Matrix_CSR>::submatrix(*matrix, A22t_dn, size_A11, size, size_A11, size, true, false, true);
		SpBLAS<T, Matrix_CSR>::submatrix(*matrix, A12t_dn, 0, size_A11, size_A11, size, true, false, true);

		cholmod_common &cm_common = _solver->cholmod.common;
		cholmod_sparse *cm_A11_sp = new cholmod_sparse();
		cholmod_sparse *cm_A21_sp = new cholmod_sparse();
		cholmod_dense *cm_A22_dn = new cholmod_dense();
		cholmod_dense *cm_A12_dn = new cholmod_dense();
		cholmod_factor *cm_L;
		cholmod_dense *cm_A11iA12_dn;

		setSymmetric(cm_A11_sp, A11_sp);
		updateSymmetric(cm_A11_sp, A11_sp);
		setSymmetric(cm_A21_sp, A21t_sp);
		updateSymmetric(cm_A21_sp, A21t_sp);
		update(cm_A22_dn, A22t_dn);
		update(cm_A12_dn, A12t_dn);

		double alpha[2] = {-1,0};
		double beta[2] = {1,0};

		_analyze<esint>(cm_L, cm_A11_sp, cm_common);
		_factorize<esint>(cm_L, cm_A11_sp, cm_common);
		_solve<esint>(cm_A11iA12_dn, cm_L, cm_A12_dn, cm_common);
		_apply<esint>(cm_A22_dn, cm_A21_sp, cm_A11iA12_dn, alpha, beta, cm_common);

		if constexpr (std::is_same_v<T,double>) { sc.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
		if constexpr (std::is_same_v<T,std::complex<double>>) { sc.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE; }
		sc.shape = Matrix_Shape::UPPER;
		sc.resize(A22t_dn);
		for(esint r = 0, i = 0; r < sc.nrows; ++r) {
			for(esint c = r; c < sc.ncols; ++c, ++i) {
				sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
			}
		}

		delete cm_A11_sp;
		delete cm_A21_sp;
		delete cm_A22_dn;
		delete cm_A12_dn;
		_free<esint>(cm_L, cm_common);
		_free<esint>(cm_A11iA12_dn, cm_common);
		break;
	}
	default:
		break; // UMFPACK
	}
}

}

#endif
#endif

