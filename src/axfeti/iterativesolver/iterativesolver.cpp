
#include "iterativesolver.h"

#include "cpg.h"
#include "pcpg.h"

namespace espreso {

template <typename T>
static IterativeSolver<T>* _set(AX_FETI<T> *feti)
{
	switch (feti->configuration.iterative_solver) {
	case FETIConfiguration::ITERATIVE_SOLVER::PCG:
		if (feti->configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE) {
			eslog::info(" = ITERATIVE SOLVER                                             CONJUGATE PROJECTED GRADIENT = \n");
			return new CPG<T>(feti);
		} else {
			eslog::info(" = ITERATIVE SOLVER                              PRECONDITIONED CONJUGATE PROJECTED GRADIENT = \n");
			return new PCPG<T>(feti);
		}
	case FETIConfiguration::ITERATIVE_SOLVER::pipePCG:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG:
	case FETIConfiguration::ITERATIVE_SOLVER::GMRES:
	case FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB:
	case FETIConfiguration::ITERATIVE_SOLVER::QPCE:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP:
	case FETIConfiguration::ITERATIVE_SOLVER::PCG_CP:
	default: return nullptr;
	}
}

template <typename T>
void _reconstructSolution(IterativeSolver<T> *solver, const Vector_Dual<T> &l, const Vector_Dual<T> &r)
{

}

template <> IterativeSolver<double>* IterativeSolver<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> IterativeSolver<std::complex<double> >* IterativeSolver<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

template <> void IterativeSolver<double>::reconstructSolution(const Vector_Dual<double> &l, const Vector_Dual<double> &r) { return _reconstructSolution<double>(this, l, r); }
template <> void IterativeSolver<std::complex<double> >::reconstructSolution(const Vector_Dual<std::complex<double> > &l, const Vector_Dual<std::complex<double> > &r) { return _reconstructSolution<std::complex<double> >(this, l, r); }

}
