
#ifndef SRC_ANALYSIS_LINEARSYSTEM_FETI_FETISYSTEM_H_
#define SRC_ANALYSIS_LINEARSYSTEM_FETI_FETISYSTEM_H_

#include "analysis/linearsystem/linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/acoustic.real.linear.h"
#include "analysis/analysis/acoustic.complex.linear.h"
#include "analysis/analysis/structuralmechanics.steadystate.linear.h"
#include "analysis/analysis/structuralmechanics.steadystate.nonlinear.h"

#include "analysis/composer/nodes.uniform.feti.h"
#include "feti/feti.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"

namespace espreso {

template <typename Assembler, typename Solver>
struct FETISystemData: public LinearSystem<Assembler, Solver> {

	FETISystemData(FETIConfiguration &configuration): feti(configuration) {}

	void setMapping(Matrix_Base<Assembler> *A) const
	{
		assembler.pattern.setMap(dynamic_cast<Matrix_FETI<Matrix_CSR, Assembler>*>(A));
	}

	void setMapping(Vector_Base<Assembler> *x) const
	{
		assembler.pattern.setMap(dynamic_cast<Vector_FETI<Vector_Dense, Assembler>*>(x));
	}

	void setDirichletMapping(Vector_Base<Assembler> *x) const
	{
		assembler.pattern.setDirichletMap(dynamic_cast<Vector_Distributed<Vector_Sparse, Assembler>*>(x));
	}

	bool solve(step::Step &step)
	{
		eslog::startln("FETI: RUN LINEAR SYSTEM", "FETI[SOLVE]");
		if (feti.solve(step, solver.x)) {
			if (info::ecf->output.print_matrices) {
				eslog::storedata(" STORE: system/{x}\n");
				math::store(solver.x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
			}
			eslog::endln("FETI: LINEAR SYSTEM SOLVED");
			return true;
		}
		eslog::endln("FETI: LINEAR SYSTEM SOLVED");
		return false;
	}

	template <typename Type>
	struct Data {
		UniformNodesFETIPattern pattern;

		Matrix_FETI<Matrix_CSR, Type> K;
		Vector_FETI<Vector_Dense, Type> x, f;
		Vector_Distributed<Vector_Sparse, Type> dirichlet;

		DOFsDecomposition decomposition;
	};

	Data<Assembler> assembler;
	Data<Solver> solver;

	FETI<Solver> feti;
	typename FETI<Solver>::Regularization regularization;
	typename FETI<Solver>::EqualityConstraints equalityConstraints;
};

template <typename T>
void composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange);
template <typename T>
void evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange);

void composeHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat);
void evaluateHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat);

template <typename Assembler, typename Solver>
void setEqualityConstraints(FETISystemData<Assembler, Solver> *system, step::Step &step)
{
	if (system->feti.configuration.redundant_lagrange) {
		eslog::info(" = LAGRANGE MULTIPLIERS                                                            REDUNDANT = \n");
	} else {
		eslog::info(" = LAGRANGE MULTIPLIERS                                                          ORTHONORMAL = \n");
	}
	composeEqualityConstraints(system->solver.K, system->solver.dirichlet, system->equalityConstraints, system->feti.configuration.redundant_lagrange);
}

template <typename Assembler, typename Solver>
void setHeatTransferKernel(FETISystemData<Assembler, Solver> *system, step::Step &step)
{
	eslog::info(" = REGULARIZATION                                                                   ANALYTIC = \n");
	system->regularization.R1.domains.resize(system->solver.K.domains.size());
	system->regularization.R2.domains.resize(system->solver.K.domains.size());
	system->regularization.RegMat.domains.resize(system->solver.K.domains.size());

	system->regularization.RegMat.type = system->solver.K.type;
	system->regularization.RegMat.shape = system->solver.K.shape;
	#pragma omp parallel for
	for (size_t d = 0; d < system->solver.K.domains.size(); ++d) {
		// TODO: some domains can be without kernel
		composeHeatTransferKernel(system->solver.K.domains[d], system->regularization.R1.domains[d], system->regularization.RegMat.domains[d]);
	}
}

template <typename Assembler, typename Solver>
void evaluateEqualityConstraints(FETISystemData<Assembler, Solver> *system, step::Step &step)
{
	evaluateEqualityConstraints(system->solver.K, system->solver.dirichlet, system->equalityConstraints, system->feti.configuration.redundant_lagrange);
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: system/{B1, B1c, B1Duplication, D2C, LMAP}\n");
		math::store(system->equalityConstraints.c, utils::filename(utils::debugDirectory(step) + "/system", "B1c").c_str());
		for (size_t d = 0; d < system->equalityConstraints.domain.size(); ++d) {
			math::store(system->equalityConstraints.domain[d].B1, utils::filename(utils::debugDirectory(step) + "/system", "B1" + std::to_string(d)).c_str());
			math::store(system->equalityConstraints.domain[d].duplication, utils::filename(utils::debugDirectory(step) + "/system", "B1Duplication" + std::to_string(d)).c_str());
			math::store(system->equalityConstraints.domain[d].D2C, utils::filename(utils::debugDirectory(step) + "/system", "D2C" + std::to_string(d)).c_str());
			math::store(system->equalityConstraints.lmap, utils::filename(utils::debugDirectory(step) + "/system", "LMAP").c_str());
		}
	}
}

template <typename Assembler, typename Solver>
void evaluateHeatTransferKernel(FETISystemData<Assembler, Solver> *system, step::Step &step)
{
	#pragma omp parallel for
	for (size_t d = 0; d < system->solver.K.domains.size(); ++d) {
		// TODO: some domains can be without kernel
		evaluateHeatTransferKernel(system->solver.K.domains[d], system->regularization.R1.domains[d], system->regularization.RegMat.domains[d]);
	}
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: system/{R, RegMat}\n");
		math::store(system->regularization.R1, utils::filename(utils::debugDirectory(step) + "/system", "R").c_str());
		math::store(system->regularization.RegMat, utils::filename(utils::debugDirectory(step) + "/system", "RegMat").c_str());
	}
}

template <typename Analysis> struct FETISystem {};

template <typename A, typename S>
inline void _fillDirect(FETISystemData<A, S> *system, std::map<std::string, ECFExpression> &dirichlet, int dofs)
{
	system->assembler.pattern.set(dirichlet, dofs, system->solver.decomposition, system->solver.K.shape);
	system->assembler.pattern.fill(system->solver.K);
	system->assembler.pattern.fill(system->solver.f);
	system->assembler.pattern.fill(system->solver.x);
	system->assembler.pattern.fill(system->solver.dirichlet);

	system->assembler.dirichlet.distribution = system->assembler.K.decomposition = system->assembler.f.decomposition = system->assembler.x.decomposition = &system->solver.decomposition;
	system->solver.dirichlet.distribution = system->solver.K.decomposition = system->solver.f.decomposition = system->solver.x.decomposition = &system->solver.decomposition;

	system->LinearSystem<A, S>::assembler.A = system->LinearSystem<A, S>::solver.A = &system->solver.K;
	system->LinearSystem<A, S>::assembler.b = system->LinearSystem<A, S>::solver.b = &system->solver.f;
	system->LinearSystem<A, S>::assembler.x = system->LinearSystem<A, S>::solver.x = &system->solver.x;
	system->LinearSystem<A, S>::assembler.dirichlet = system->LinearSystem<A, S>::solver.dirichlet = &system->solver.dirichlet;
}

template <> struct FETISystem<HeatSteadyStateLinear>: public FETISystemData<double, double> {

	FETISystem(HeatSteadyStateLinear *analysis): FETISystemData(analysis->configuration.feti)
	{
		assembler.K.type = solver.K.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
		assembler.K.shape = solver.K.shape = Matrix_Shape::UPPER;
		for (auto mat = analysis->settings.material_set.begin(); mat != analysis->settings.material_set.end(); ++mat) {
			if (analysis->settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
				assembler.K.type = solver.K.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
				assembler.K.shape = solver.K.shape = Matrix_Shape::FULL;
			}
		}

		_fillDirect(this, analysis->configuration.temperature, 1);
	}

	void set(step::Step &step)
	{
		eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
		setEqualityConstraints(this, step);
		eslog::checkpointln("FETI: SET B1");
		setHeatTransferKernel(this, step);
		eslog::checkpointln("FETI: SET KERNELS");
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
		feti.set(step, solver.K, regularization, equalityConstraints);
		eslog::endln("FETI: LINEAR SYSTEM SET");
	}

	void update(step::Step &step)
	{
		eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
		evaluateEqualityConstraints(this, step);
		eslog::checkpointln("FETI: UPDATE B1");
		evaluateHeatTransferKernel(this, step);
		eslog::checkpointln("FETI: UPDATE KERNELS");
		feti.update(step, solver.K, solver.f);
		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: system/{K, f}\n");
			math::store(solver.K, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
			math::store(solver.f, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
		}
		eslog::endln("FETI: LINEAR SYSTEM UPDATED");
	}
};

template <> struct FETISystem<HeatSteadyStateNonLinear>: public FETISystemData<double, double> {

	FETISystem(HeatSteadyStateNonLinear *analysis): FETISystemData(analysis->configuration.feti)
	{

	}

	void set(step::Step &step)
	{
		setEqualityConstraints(this, step);
	}

	void update(step::Step &step)
	{

	}
};

template <> struct FETISystem<StructuralMechanicsSteadyStateLinear>: public FETISystemData<double, double> {

	FETISystem(StructuralMechanicsSteadyStateLinear *analysis): FETISystemData(analysis->configuration.feti)
	{
//		assembler.K.type = solver.K.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//		assembler.K.shape = solver.K.shape = Matrix_Shape::UPPER;
//		for (auto mat = analysis->settings.material_set.begin(); mat != analysis->settings.material_set.end(); ++mat) {
//			if (analysis->settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
//				assembler.K.type = solver.K.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
//				assembler.K.shape = solver.K.shape = Matrix_Shape::FULL;
//			}
//		}
//
//		_fillDirect(this, analysis->configuration.temperature, 1);
	}

	void set(step::Step &step)
	{
//		eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
//		setEqualityConstraints(this, step);
//		eslog::checkpointln("FETI: SET B1");
//		setHeatTransferKernel(this, step);
//		eslog::checkpointln("FETI: SET KERNELS");
//		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
//		feti.set(step, solver.K, regularization, equalityConstraints);
//		eslog::endln("FETI: LINEAR SYSTEM SET");
	}

	void update(step::Step &step)
	{
//		eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
//		evaluateEqualityConstraints(this, step);
//		eslog::checkpointln("FETI: UPDATE B1");
//		evaluateHeatTransferKernel(this, step);
//		eslog::checkpointln("FETI: UPDATE KERNELS");
//		feti.update(step, solver.K, solver.f);
//		if (info::ecf->output.print_matrices) {
//			eslog::storedata(" STORE: system/{K, f}\n");
//			math::store(solver.K, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
//			math::store(solver.f, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
//		}
//		eslog::endln("FETI: LINEAR SYSTEM UPDATED");
	}
};

template <> struct FETISystem<StructuralMechanicsSteadyStateNonLinear>: public FETISystemData<double, double> {

	FETISystem(StructuralMechanicsSteadyStateNonLinear *analysis): FETISystemData(analysis->configuration.feti)
	{
//		assembler.K.type = solver.K.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//		assembler.K.shape = solver.K.shape = Matrix_Shape::UPPER;
//		for (auto mat = analysis->settings.material_set.begin(); mat != analysis->settings.material_set.end(); ++mat) {
//			if (analysis->settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
//				assembler.K.type = solver.K.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
//				assembler.K.shape = solver.K.shape = Matrix_Shape::FULL;
//			}
//		}
//
//		_fillDirect(this, analysis->configuration.temperature, 1);
	}

	void set(step::Step &step)
	{
//		eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
//		setEqualityConstraints(this, step);
//		eslog::checkpointln("FETI: SET B1");
//		setHeatTransferKernel(this, step);
//		eslog::checkpointln("FETI: SET KERNELS");
//		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
//		feti.set(step, solver.K, regularization, equalityConstraints);
//		eslog::endln("FETI: LINEAR SYSTEM SET");
	}

	void update(step::Step &step)
	{
//		eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
//		evaluateEqualityConstraints(this, step);
//		eslog::checkpointln("FETI: UPDATE B1");
//		evaluateHeatTransferKernel(this, step);
//		eslog::checkpointln("FETI: UPDATE KERNELS");
//		feti.update(step, solver.K, solver.f);
//		if (info::ecf->output.print_matrices) {
//			eslog::storedata(" STORE: system/{K, f}\n");
//			math::store(solver.K, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
//			math::store(solver.f, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
//		}
//		eslog::endln("FETI: LINEAR SYSTEM UPDATED");
	}
};

template <> struct FETISystem<AcousticRealLinear>: public FETISystemData<double, double> {

	FETISystem(AcousticRealLinear *analysis): FETISystemData(analysis->configuration.feti)
	{

	}

	void set(step::Step &step)
	{
		setEqualityConstraints(this, step);
	}

	void update(step::Step &step)
	{

	}
};

template <> struct FETISystem<AcousticComplexLinear>: public FETISystemData<double, std::complex<double> > {

	FETISystem(AcousticComplexLinear *analysis)
	: FETISystemData(analysis->configuration.feti)
	{

	}

	void set(step::Step &step)
	{
		setEqualityConstraints(this, step);
	}

	void update(step::Step &step)
	{

	}
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_FETI_FETISYSTEM_H_ */
