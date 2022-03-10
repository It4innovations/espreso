
#include "iterativesolver.h"

#include "pcpg.h"

namespace espreso {

template <typename T>
static IterativeSolver<T>* _set(AX_FETI<T> *feti)
{
	switch (feti->configuration.iterative_solver) {
	case FETIConfiguration::ITERATIVE_SOLVER::PCG:
		eslog::info(" = ITERATIVE SOLVER                                                                      PCG = \n");
		return new PCPG<T>(feti);
		break;
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

template <> IterativeSolver<double>* IterativeSolver<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> IterativeSolver<std::complex<double> >* IterativeSolver<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
