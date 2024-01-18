
#include "fetisolver.h"

#include "regularization/regularization.heattransfer.h"
#include "regularization/regularization.elasticity.h"
#include "constrains/equalityconstrains.h"

#include "analysis/physics/heat.steadystate.linear.h"
#include "analysis/physics/heat.steadystate.nonlinear.h"
#include "analysis/physics/structuralmechanics.steadystate.linear.h"
#include "analysis/physics/structuralmechanics.steadystate.nonlinear.h"

namespace espreso {

template <typename T, class Physics> struct RegularizationSelector { };

template <typename T> struct RegularizationSelector<T, HeatSteadyStateLinear> {
	Regularization<T>* operator()(FETI<T> &feti){ return new RegularizationHeatTransfer<T>(feti); }
};

template <typename T> struct RegularizationSelector<T, HeatSteadyStateNonLinear> {
	Regularization<T>* operator()(FETI<T> &feti){ return new RegularizationHeatTransfer<T>(feti); }
};

template <typename T> struct RegularizationSelector<T, StructuralMechanicsSteadyStateLinear> {
	Regularization<T>* operator()(FETI<T> &feti){ return new RegularizationElasticity<T>(feti); }
};

template <typename T> struct RegularizationSelector<T, StructuralMechanicsSteadyStateNonLinear> {
	Regularization<T>* operator()(FETI<T> &feti){ return new RegularizationElasticity<T>(feti); }
};

template <typename T, class Physics>
FETILinearSystemSolver<T, Physics>::FETILinearSystemSolver(FETIConfiguration &configuration)
: feti(configuration),
  equalityConstrains(new EqualityConstrains<T>(feti)),
  regularization(RegularizationSelector<T, Physics>()(feti))
{
	LinearSystemSolver<T>::A = &feti.K;
	LinearSystemSolver<T>::x = &feti.x;
	LinearSystemSolver<T>::b = &feti.f;
	LinearSystemSolver<T>::dirichlet = &dirichlet;
}

template <typename T, class Physics>
FETILinearSystemSolver<T, Physics>::~FETILinearSystemSolver()
{
	delete equalityConstrains;
	delete regularization;
}

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::set(step::Step &step)
{
	eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
	feti.decomposition = feti.K.decomposition;
	equalityConstrains->set(step, dirichlet);
	eslog::checkpointln("FETI: SET B1");
	regularization->set(step);
	eslog::checkpointln("FETI: SET KERNELS");
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	feti.set(step);
	eslog::endln("FETI: LINEAR SYSTEM SET");
}

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::update(step::Step &step)
{
	eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
	equalityConstrains->update(step, dirichlet);
	eslog::checkpointln("FETI: UPDATE B1");
	regularization->update(step);
	eslog::checkpointln("FETI: UPDATE KERNELS");
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: system/{K, f, R, RegMat}\n");
		math::store(feti.K, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
		math::store(feti.f, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
		for (size_t d = 0; d < feti.regularization.R1.size(); ++d) {
			math::store(feti.regularization.R1[d], utils::filename(utils::debugDirectory(step) + "/system", "R" + std::to_string(d)).c_str());
			math::store(feti.regularization.RegMat[d], utils::filename(utils::debugDirectory(step) + "/system", "RegMat" + std::to_string(d)).c_str());
		}

		eslog::storedata(" STORE: system/{B1, B1c, B1Duplication, D2C, LMAP}\n");
		math::store(feti.equalityConstraints.c, utils::filename(utils::debugDirectory(step) + "/system", "B1c").c_str());
		for (size_t d = 0; d < feti.equalityConstraints.domain.size(); ++d) {
			math::store(feti.equalityConstraints.domain[d].B1, utils::filename(utils::debugDirectory(step) + "/system", "B1" + std::to_string(d)).c_str());
			math::store(feti.equalityConstraints.domain[d].D2C, utils::filename(utils::debugDirectory(step) + "/system", "D2C" + std::to_string(d)).c_str());
		}
	}
	feti.update(step);
	eslog::endln("FETI: LINEAR SYSTEM UPDATED");
}

template <typename T, class Physics>
bool FETILinearSystemSolver<T, Physics>::solve(step::Step &step)
{
	eslog::startln("FETI: RUN LINEAR SYSTEM", "FETI[SOLVE]");
	if (feti.solve(step)) {
		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: system/{x}\n");
			math::store(feti.x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
		}
		eslog::endln("FETI: LINEAR SYSTEM SOLVED");
		return true;
	}
	eslog::endln("FETI: LINEAR SYSTEM SOLVED");
	return false;
}

template struct FETILinearSystemSolver<double, HeatSteadyStateLinear>;
template struct FETILinearSystemSolver<double, HeatSteadyStateNonLinear>;
template struct FETILinearSystemSolver<double, StructuralMechanicsSteadyStateLinear>;
template struct FETILinearSystemSolver<double, StructuralMechanicsSteadyStateNonLinear>;

}

