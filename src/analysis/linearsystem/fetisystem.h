
#ifndef SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/acoustic.real.linear.h"
#include "analysis/analysis/acoustic.complex.linear.h"

#include "analysis/composer/nodes.uniform.feti.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "axfeti/feti.h"

namespace espreso {

template <typename Assembler, typename Solver>
struct AX_FETISystemData: public AX_LinearSystem<Assembler, Solver> {

	AX_FETISystemData(FETIConfiguration &configuration) {}

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

	void info() const
	{
//		feti.info(solver.A);
	}

	void set(step::Step &step)
	{
//		feti.set(solver.A);
	}

	void update(step::Step &step)
	{

	}

	bool solve(step::Step &step)
	{

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

	AX_FETI<Solver> feti;
};

template <typename Analysis> struct AX_FETISystem {};

template <typename A, typename S>
inline void _fillDirect(AX_FETISystemData<A, S> *system, std::map<std::string, ECFExpression> &dirichlet, int dofs)
{
	system->assembler.pattern.set(dirichlet, dofs, system->solver.decomposition, system->solver.K.shape);
	system->assembler.pattern.fill(system->solver.K);
	system->assembler.pattern.fill(system->solver.f);
	system->assembler.pattern.fill(system->solver.x);
	system->assembler.pattern.fill(system->solver.dirichlet);

	system->assembler.dirichlet.distribution = system->assembler.K.decomposition = system->assembler.f.decomposition = system->assembler.x.decomposition = &system->solver.decomposition;
	system->solver.dirichlet.distribution = system->solver.K.decomposition = system->solver.f.decomposition = system->solver.x.decomposition = &system->solver.decomposition;

	system->AX_LinearSystem<A, S>::assembler.A = system->AX_LinearSystem<A, S>::solver.A = &system->solver.K;
	system->AX_LinearSystem<A, S>::assembler.b = system->AX_LinearSystem<A, S>::solver.b = &system->solver.f;
	system->AX_LinearSystem<A, S>::assembler.x = system->AX_LinearSystem<A, S>::solver.x = &system->solver.x;
	system->AX_LinearSystem<A, S>::assembler.dirichlet = system->AX_LinearSystem<A, S>::solver.dirichlet = &system->solver.dirichlet;
}

template <> struct AX_FETISystem<AX_HeatSteadyStateLinear>: public AX_FETISystemData<double, double> {

	AX_FETISystem(AX_HeatSteadyStateLinear *analysis)
	: AX_FETISystemData(analysis->configuration.feti)
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
};

template <> struct AX_FETISystem<AX_HeatSteadyStateNonLinear>: public AX_FETISystemData<double, double> {

	AX_FETISystem(AX_HeatSteadyStateNonLinear *analysis)
	: AX_FETISystemData(analysis->configuration.feti)
	{

	}
};

template <> struct AX_FETISystem<AX_AcousticRealLinear>: public AX_FETISystemData<double, double> {

	AX_FETISystem(AX_AcousticRealLinear *analysis)
	: AX_FETISystemData(analysis->configuration.feti)
	{

	}
};

template <> struct AX_FETISystem<AX_AcousticComplexLinear>: public AX_FETISystemData<double, std::complex<double> > {

	AX_FETISystem(AX_AcousticComplexLinear *analysis)
	: AX_FETISystemData(analysis->configuration.feti)
	{

	}
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_ */
