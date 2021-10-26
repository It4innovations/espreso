
#include "feti.h"
#include "iterativesolver/pcg.h"
#include "preconditioner/emptypreconditioner.h"

namespace espreso {

template <typename T>
bool _set(AX_FETI<T> *feti, const Matrix_FETI<Matrix_CSR, T> &K, const typename AX_FETI<T>::Regularization &regularization, const typename AX_FETI<T>::EqualityConstraints &equalityConstraints)
{
	switch (feti->configuration.iterative_solver) {
	case FETIConfiguration::ITERATIVE_SOLVER::PCG: feti->iterativeSolver = new PCG<T>(); break;
	case FETIConfiguration::ITERATIVE_SOLVER::pipePCG:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG:
	case FETIConfiguration::ITERATIVE_SOLVER::GMRES:
	case FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB:
	case FETIConfiguration::ITERATIVE_SOLVER::QPCE:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP:
	case FETIConfiguration::ITERATIVE_SOLVER::PCG_CP:
	default:
		feti->iterativeSolver = new IterativeSolver<T>();
	}

	switch (feti->configuration.preconditioner) {
	case FETIConfiguration::PRECONDITIONER::NONE: feti->preconditioner = new EmptyPreconditioner<T>(); break;
	case FETIConfiguration::PRECONDITIONER::LUMPED:
	case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
	case FETIConfiguration::PRECONDITIONER::DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::MAGIC:
	default:
		feti->preconditioner = new Preconditioner<T>();
	}

//	switch (feti->configuration.projector) {
//
//	}

	// G = R'B' (do we need to know values in R?)

	// 1. per process R'B' -> rows of G
	// 2. Gather G
	// 3. compute GG'
	// 4a. Factorization of GG'
	// 4b. Explicit inv(GG')
	return true;
}

template <typename T>
bool _update(const Matrix_Distributed<Matrix_CSR, T> &K)
{
	return true;
}

template <typename T>
bool _solve(const Vector_Distributed<Vector_Dense, T> &f, Vector_Distributed<Vector_Dense, T> &x)
{
	return true;
}


template <> bool AX_FETI<double>::set(const Matrix_FETI<Matrix_CSR, double> &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints) { return _set(this, K, regularization, equalityConstraints); }
template <> bool AX_FETI<std::complex<double> >::set(const Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints) { return _set(this, K, regularization, equalityConstraints); }

template <> bool AX_FETI<double>::update(const Matrix_Distributed<Matrix_CSR, double> &K) { return _update(K); }
template <> bool AX_FETI<std::complex<double> >::update(const Matrix_Distributed<Matrix_CSR, std::complex<double> > &K) { return _update(K); }

template <> bool AX_FETI<double>::solve(const Vector_Distributed<Vector_Dense, double> &f, Vector_Distributed<Vector_Dense, double> &x) { return _solve(f, x); }
template <> bool AX_FETI<std::complex<double> >::solve(const Vector_Distributed<Vector_Dense, std::complex<double> > &f, Vector_Distributed<Vector_Dense, std::complex<double> > &x) { return _solve(f, x); }

}
