
#include "math/math.h"
#include "esinfo/eslog.h"

#ifdef HAVE_MKL
#ifdef USE_SOLVER_MKL

#include "wrappers/pardiso/w.pardiso.type.h"
#include "wrappers/pardiso/w.pardiso.h"
#include "mkl_pardiso.h"

namespace espreso {

struct Matrix_CSR_Solver: public PARDISOParameters { };

namespace math {

template<typename T>
bool _callPardiso(esint phase, const Matrix_CSR<T> &m, esint nrhs, T *rhs, T *solution)
{
	m._solver->phase = phase;
	pardiso(
			m._solver->pt, &m._solver->maxfct, &m._solver->mnum,
			&m._solver->mtype,
			&m._solver->phase,
			&m.nrows, m.vals, m.rows, m.cols,
			m._solver->perm, &nrhs, m._solver->iparm, &m._solver->msglvl,
			rhs, solution,
			&m._solver->error);

	switch (m._solver->error) {
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
	return m._solver->error == 0;
}

template <>
void initSolver(Matrix_CSR<double> &x)
{
	x._solver = new Matrix_CSR_Solver();
	x._solver->mtype = _pardisoType(x);
	pardisoinit(x._solver->pt, &x._solver->mtype, x._solver->iparm);
	x._solver->iparm[0] = 1;			/* No solver default */
	x._solver->iparm[1] = 2;			/* Fill-in reordering from METIS */
	x._solver->iparm[2] = 1; 			//by default the solver runs with single thread
	x._solver->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	x._solver->iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
}

template <>
void initSolver(Matrix_CSR<std::complex<double> > &x)
{
	x._solver = new Matrix_CSR_Solver();
	x._solver->mtype = _pardisoType(x);
	pardisoinit(x._solver->pt, &x._solver->mtype, x._solver->iparm);
	x._solver->iparm[0] = 1;			/* No solver default */
	x._solver->iparm[1] = 2;			/* Fill-in reordering from METIS */
	x._solver->iparm[2] = 1; 			//by default the solver runs with single thread
	x._solver->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	x._solver->iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
}

template <>
void symbolicFactorization(const Matrix_CSR<double> &x)
{
	_callPardiso<double>(11, x, 0, nullptr, nullptr);
}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &x)
{
	x._solver->mtype = _pardisoType(x);
	pardisoinit(x._solver->pt, &x._solver->mtype, x._solver->iparm);
	x._solver->iparm[0] = 1;			/* No solver default */
	x._solver->iparm[1] = 2;			/* Fill-in reordering from METIS */
	x._solver->iparm[2] = 1; 			//by default the solver runs with single thread
	x._solver->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	x._solver->iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
	_callPardiso<std::complex<double> >(11, x, 0, nullptr, nullptr);
}

template <>
void numericalFactorization(const Matrix_CSR<double> &x)
{
	_callPardiso<double>(22, x, 0, nullptr, nullptr);
}

template <>
void numericalFactorization(const Matrix_CSR<std::complex<double> > &x)
{
	_callPardiso<std::complex<double> >(22, x, 0, nullptr, nullptr);
}

template <>
void solve(const Matrix_CSR<double> &x, Vector_Dense<double> &rhs, Vector_Dense<double> &solution)
{
	_callPardiso<double>(33, x, 1, rhs.vals, solution.vals);
}

template <>
void solve(const Matrix_CSR<double> &x, Matrix_Dense<double> &rhs, Matrix_Dense<double> &solution)
{
	_callPardiso<double>(33, x, rhs.nrows, rhs.vals, solution.vals);
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Vector_Dense<std::complex<double> > &rhs, Vector_Dense<std::complex<double> > &solution)
{
	_callPardiso<std::complex<double> >(33, x, 1, rhs.vals, solution.vals);
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Matrix_Dense<std::complex<double> > &rhs, Matrix_Dense<std::complex<double> > &solution)
{
	_callPardiso<std::complex<double> >(33, x, rhs.nrows, rhs.vals, solution.vals);
}

template <>
void freeSolver(Matrix_CSR<double> &x)
{
	_callPardiso<double>(-1, x, 0, nullptr, nullptr);
	delete x._solver;
}

template <>
void freeSolver(Matrix_CSR<std::complex<double> > &x)
{
	_callPardiso<std::complex<double> >(-1, x, 0, nullptr, nullptr);
	delete x._solver;
}

}
}

#endif
#endif // HAVE_MKL
