
#include "smalbe.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

namespace espreso {

template <typename T>
SMALBE<T>::SMALBE(FETI<T> &feti)
: IterativeSolver<T>(feti)
{

}

template <typename T>
void SMALBE<T>::info()
{
	eslog::info(" = CONJUGATE PROJECTED GRADIENT SETTINGS                                                     = \n");
	switch (feti.configuration.stopping_criterion) {
	case FETIConfiguration::STOPPING_CRITERION::RELATIVE:
		eslog::info(" =   STOPPING CRITERION                                                             RELATIVE = \n");
		break;
	case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE:
		eslog::info(" =   STOPPING CRITERION                                                             ABSOLUTE = \n");
		break;
	case FETIConfiguration::STOPPING_CRITERION::ARIOLI:
		eslog::info(" =   STOPPING CRITERION                                                               ARIOLI = \n");
		break;
	}
	eslog::info(" =   PRECISION                                                                      %.2e = \n", feti.configuration.precision);
	if (feti.configuration.max_iterations == 0) {
		eslog::info(" =   MAX_ITERATIONS                                                                     AUTO = \n");
	} else {
		eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", feti.configuration.max_iterations);
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _print(const char *name, const IterativeSolverInfo &info, const step::Step &step, const Vector_Dual<T> &v)
{
	if (info::ecf->output.print_matrices > 1) {
		eslog::storedata(" STORE: feti/SMALBE/{%s%s}\n", name, std::to_string(info.iterations).c_str());
		math::store(v, utils::filename(utils::debugDirectory(step) + "/feti/SMALBE", std::string(name) + std::to_string(info.iterations)).c_str());
	}
}

template <> void SMALBE<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template <> void SMALBE<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class SMALBE<double>;
template class SMALBE<std::complex<double> >;

}