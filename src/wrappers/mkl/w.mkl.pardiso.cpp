
#include "math/math.h"
#include "esinfo/eslog.h"

#ifdef HAVE_MKL
#ifdef USE_SOLVER_MKL

#include "wrappers/pardiso/w.pardiso.type.h"
#include "wrappers/pardiso/w.pardiso.h"
#include "mkl_pardiso.h"

namespace espreso {

template class DirectSolver<double, Matrix_CSR>;
template class DirectSolver<std::complex<double>, Matrix_CSR>;

struct Matrix_CSR_Solver: public PARDISOParameters { };

struct Matrix_Solver_External_Representation: PARDISOParameters { };

template<typename T>
bool _callPardiso(esint phase, const Matrix_CSR<T> &m, Matrix_Solver_External_Representation *solver, esint nrhs, T *rhs, T *solution)
{
	solver->phase = phase;
	pardiso(
			solver->pt, &solver->maxfct, &solver->mnum,
			&solver->mtype,
			&solver->phase,
			&m.nrows, m.vals, m.rows, m.cols,
			solver->perm, &nrhs, solver->iparm, &solver->msglvl,
			rhs, solution,
			&solver->error);

	switch (solver->error) {
	case   0: break;
	case  -1: eslog::error("MKL PARDISO: input inconsistent.\n"); break;
	case  -2: eslog::error("MKL PARDISO: not enough memory.\n"); break;
	case  -3: eslog::error("MKL PARDISO: reordering problem.\n"); break;
	case  -4: eslog::error("MKL PARDISO: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
	case  -5: eslog::error("MKL PARDISO: unclassified (internal) error.\n"); break;
	case  -6: eslog::error("MKL PARDISO: reordering failed.\n"); break;
	case  -7: eslog::error("MKL PARDISO: diagonal matrix is singular.\n"); break;
	case  -8: eslog::error("MKL PARDISO: 32-bit integer overflow problem.\n"); break;
	case  -9: eslog::error("MKL PARDISO: not enough memory for OOC.\n"); break;
	case -10: eslog::error("MKL PARDISO: error opening OOC files.\n"); break;
	case -11: eslog::error("MKL PARDISO: read/write error with OOC files.\n"); break;
	case -12: eslog::error("MKL PARDISO: (pardiso_64 only) pardiso_64 called from 32-bit library.\n"); break;
	case -13: eslog::error("MKL PARDISO: interrupted by the (user-defined) mkl_progress function.\n"); break;
	case -15: eslog::error("MKL PARDISO: internal error which can appear for iparm[23]=10 and iparm[12]=1. Try switch matching off (set iparm[12]=0 and rerun.\n"); break;
	}
	return solver->error == 0;
}

template <typename T, template <typename> class Matrix>
const char* DirectSolver<T, Matrix>::name()
{
	return "MKL PARDISO";
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideFactors()
{
	return false;
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideSC()
{
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
		_callPardiso<T>(-1, *matrix, _solver, 0, nullptr, nullptr);
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
	_solver = new Matrix_Solver_External_Representation();
	_solver->mtype = _pardisoType(*matrix);
	pardisoinit(_solver->pt, &_solver->mtype, _solver->iparm);
	_solver->iparm[0] = 1;			/* No solver default */
	_solver->iparm[1] = 2;			/* Fill-in reordering from METIS */
	_solver->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::symbolicFactorization(int fixedSuffix)
{
	if (fixedSuffix) {
		_solver->iparm[30] = 1;
		_solver->perm = new esint[matrix->nrows];
		for (esint i = 0          ; i < fixedSuffix;   ++i) { _solver->perm[i] = 0; }
		for (esint i = fixedSuffix; i < matrix->nrows; ++i) { _solver->perm[i] = 1; }
	}

	_callPardiso<T>(11, *matrix, _solver, 0, nullptr, nullptr);
	rows = matrix->nrows;
	nnzA = matrix->nnz;
	memoryL = 1024 * (_solver->iparm[15] + _solver->iparm[16]);
	nnzL = _solver->iparm[17];
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::numericalFactorization()
{
	_callPardiso<T>(22, *matrix, _solver, 0, nullptr, nullptr);
}

// https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface/pardiso-iparm-parameter.html
template <typename T, template <typename> class Matrix>
inline int _sparsity(int sparsity)
{
	if (sparsity == (DirectSolver<T, Matrix>::VectorSparsity::SPARSE_RHS | DirectSolver<T, Matrix>::VectorSparsity::SPARSE_SOLUTION)) {
		return 1;
	}
	if (sparsity == DirectSolver<T, Matrix>::VectorSparsity::SPARSE_RHS) {
		return 2;
	}
	if (sparsity == DirectSolver<T, Matrix>::VectorSparsity::SPARSE_RHS) {
		return 3;
	}
	return 0;
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Vector_Dense<T> &rhs, Vector_Dense<T> &solution, int sparsity)
{
	_solver->iparm[30] = _sparsity<T, Matrix>(sparsity);
	_callPardiso<T>(33, *matrix, _solver, 1, rhs.vals, solution.vals);
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution, int sparsity)
{
	_solver->iparm[30] = _sparsity<T, Matrix>(sparsity);
	_callPardiso<T>(33, *matrix, _solver, rhs.nrows, rhs.vals, solution.vals);
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::getFactors(Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p)
{
	eslog::error("MKL PARDISO does not provide factors.\n");
}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::getSC(Matrix_Dense<T> &sc)
{
	Matrix_Dense<T> full; full.resize(sc.nrows, sc.nrows);
	Vector_Dense<esint> perm; perm.resize(matrix->nrows);
	for (esint i = 0                       ; i < matrix->nrows - sc.nrows; ++i) { perm.vals[i] = 0; }
	for (esint i = matrix->nrows - sc.nrows; i < matrix->nrows           ; ++i) { perm.vals[i] = 1; }

	_solver->iparm[35] = 1;
	std::swap(_solver->perm, perm.vals);
	_callPardiso<T>(12, *matrix, _solver, 0, nullptr, full.vals);
	std::swap(_solver->perm, perm.vals);
	_solver->iparm[35] = 0;

	for (esint r = 0, i = 0; r < full.nrows; ++r) {
		for (esint c = r; c < full.ncols; ++c, ++i) {
			sc.vals[i] = full.vals[r * full.ncols + c];
		}
	}
}

}

#endif
#endif // HAVE_MKL
