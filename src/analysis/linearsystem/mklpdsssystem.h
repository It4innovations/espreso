
#ifndef SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/acoustic.real.linear.h"
#include "analysis/analysis/acoustic.complex.linear.h"

#include "analysis/composer/nodes.uniform.distributed.h"
#include "basis/utilities/sysutils.h"
#include "config/ecf/linearsolver/mklpdss.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/utils/distributed/distribution.h"
#include "math2/utils/distributed/synchronization.h"
#include "wrappers/mklpdss/w.mkl.pdss.h"

namespace espreso {

template <typename Assembler, typename Solver>
struct AX_MKLPDSSSystemData: public AX_LinearSystem<Assembler, Solver> {

	AX_MKLPDSSSystemData(MKLPDSSConfiguration &configuration): mklpdss(configuration) {}

	void setMapping(Matrix_Base<Assembler> *A) const
	{
		assembler.pattern.setMap(dynamic_cast<Matrix_Distributed<Matrix_CSR, Assembler>*>(A));
	}

	void setMapping(Vector_Base<Assembler> *x) const
	{
		assembler.pattern.setMap(dynamic_cast<Vector_Distributed<Vector_Dense, Assembler>*>(x));
	}

	void setDirichletMapping(Vector_Base<Assembler> *x) const
	{
		assembler.pattern.setDirichletMap(dynamic_cast<Vector_Distributed<Vector_Sparse, Assembler>*>(x));
	}

	void info() const
	{
		mklpdss.info(solver.A);
	}

	void set(step::Step &step)
	{
		mklpdss.set(solver.A);
	}

	void update(step::Step &step)
	{
		if (solver.A.touched || solver.b.touched || solver.dirichlet.touched) {
			setDirichlet(solver.A, solver.b, solver.dirichlet.cluster, *solver.A.distribution);
			mklpdss.update(solver.A);
		}

		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: system/{A, b, dirichlet}\n");
			math::store(solver.A, utils::filename(utils::debugDirectory(step) + "/system", "A").c_str());
			math::store(solver.b, utils::filename(utils::debugDirectory(step) + "/system", "b").c_str());
			math::store(solver.dirichlet, utils::filename(utils::debugDirectory(step) + "/system", "dirichlet").c_str());
		}
	}

	bool solve(step::Step &step)
	{
		if (mklpdss.solve(solver.b, solver.x)) {
			solver.x.scatter();
			if (info::ecf->output.print_matrices) {
				eslog::storedata(" STORE: system/{x}\n");
				math::store(solver.x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
			}
			return true;
		}
		return false;
	}

	template <typename Type>
	struct Data {
		UniformNodesDistributedPattern pattern;

		Matrix_Distributed<Matrix_CSR, Type> A;
		Vector_Distributed<Vector_Dense, Type> x, b;
		Vector_Distributed<Vector_Sparse, Type> dirichlet;

		DOFsDistribution distribution;
		struct {
			Data_Synchronization<Matrix_CSR, Type> A;
			Data_Synchronization<Vector_Dense, Type> b;
		} sync;
	};
	
	Data<Assembler> assembler;
	Data<Solver> solver;

	MKLPDSS<Solver> mklpdss;
};

template <typename Analysis> struct AX_MKLPDSSSystem {};

template <typename A, typename S>
inline void _fillAssembler(AX_MKLPDSSSystemData<A, S> *system, std::map<std::string, ECFExpression> &dirichlet, int dofs)
{
	system->assembler.pattern.set(dirichlet, dofs, system->assembler.distribution);
	system->assembler.pattern.fill(system->assembler.A);
	system->assembler.pattern.fill(system->assembler.b);
	system->assembler.pattern.fill(system->assembler.x);
	system->assembler.pattern.fill(system->assembler.dirichlet);

	system->assembler.A.distribution = system->assembler.b.distribution = system->assembler.x.distribution = system->assembler.dirichlet.distribution = &system->assembler.distribution;
	system->assembler.sync.A.init(system->assembler.A);
	system->assembler.sync.b.init(system->assembler.b);
	system->assembler.A.synchronization = &system->assembler.sync.A;
	system->assembler.b.synchronization = system->assembler.x.synchronization = &system->assembler.sync.b;

	system->AX_LinearSystem<A, S>::assembler.A = &system->assembler.A;
	system->AX_LinearSystem<A, S>::assembler.b = &system->assembler.b;
	system->AX_LinearSystem<A, S>::assembler.x = &system->assembler.x;
	system->AX_LinearSystem<A, S>::assembler.dirichlet = &system->assembler.dirichlet;
}

template <typename A, typename S>
inline void _fillSolver(AX_MKLPDSSSystemData<A, S> *system, std::map<std::string, ECFExpression> &dirichlet, int dofs)
{
	system->solver.pattern.set(dirichlet, dofs, system->solver.distribution);
	system->solver.pattern.fill(system->solver.A);
	system->solver.pattern.fill(system->solver.b);
	system->solver.pattern.fill(system->solver.x);
	system->solver.pattern.fill(system->solver.dirichlet);

	system->solver.A.distribution = system->solver.b.distribution = system->solver.x.distribution = system->solver.dirichlet.distribution = &system->solver.distribution;
	system->solver.sync.A.init(system->solver.A);
	system->solver.sync.b.init(system->solver.b);
	system->solver.A.synchronization = &system->solver.sync.A;
	system->solver.b.synchronization = system->solver.x.synchronization = &system->solver.sync.b;

	system->AX_LinearSystem<A, S>::solver.A = &system->solver.A;
	system->AX_LinearSystem<A, S>::solver.b = &system->solver.b;
	system->AX_LinearSystem<A, S>::solver.x = &system->solver.x;
	system->AX_LinearSystem<A, S>::solver.dirichlet = &system->solver.dirichlet;
}

template <typename A, typename S>
inline void _fillDirect(AX_MKLPDSSSystemData<A, S> *system, std::map<std::string, ECFExpression> &dirichlet, int dofs)
{
	system->assembler.pattern.set(dirichlet, dofs, system->solver.distribution);
	system->assembler.pattern.fill(system->solver.A);
	system->assembler.pattern.fill(system->solver.b);
	system->assembler.pattern.fill(system->solver.x);
	system->assembler.pattern.fill(system->solver.dirichlet);

	system->assembler.A.distribution = system->assembler.b.distribution = system->assembler.x.distribution = system->assembler.dirichlet.distribution = &system->solver.distribution;
	system->solver.A.distribution = system->solver.b.distribution = system->solver.x.distribution = system->solver.dirichlet.distribution = &system->solver.distribution;
	system->solver.sync.A.init(system->solver.A);
	system->solver.sync.b.init(system->solver.b);
	system->solver.A.synchronization = &system->solver.sync.A;
	system->solver.b.synchronization = system->solver.x.synchronization = &system->solver.sync.b;
	system->assembler.A.synchronization = &system->solver.sync.A;
	system->assembler.b.synchronization = system->assembler.x.synchronization = &system->solver.sync.b;

	system->AX_LinearSystem<A, S>::assembler.A = system->AX_LinearSystem<A, S>::solver.A = &system->solver.A;
	system->AX_LinearSystem<A, S>::assembler.b = system->AX_LinearSystem<A, S>::solver.b = &system->solver.b;
	system->AX_LinearSystem<A, S>::assembler.x = system->AX_LinearSystem<A, S>::solver.x = &system->solver.x;
	system->AX_LinearSystem<A, S>::assembler.dirichlet = system->AX_LinearSystem<A, S>::solver.dirichlet = &system->solver.dirichlet;
}

template <> struct AX_MKLPDSSSystem<AX_HeatSteadyStateLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_HeatSteadyStateLinear *analysis)
	: AX_MKLPDSSSystemData(analysis->configuration.mklpdss)
	{
		assembler.A.type = solver.A.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
		for (auto mat = analysis->settings.material_set.begin(); mat != analysis->settings.material_set.end(); ++mat) {
			if (analysis->settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
				assembler.A.type = solver.A.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
			}
		}

		_fillDirect(this, analysis->configuration.temperature, 1);
	}
};

template <> struct AX_MKLPDSSSystem<AX_HeatSteadyStateNonLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_HeatSteadyStateNonLinear *analysis)
	: AX_MKLPDSSSystemData(analysis->configuration.mklpdss)
	{
		assembler.A.type = solver.A.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
		for (auto mat = analysis->settings.material_set.begin(); mat != analysis->settings.material_set.end(); ++mat) {
			if (analysis->settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
				assembler.A.type = solver.A.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
			}
		}

		_fillDirect(this, analysis->configuration.temperature, 1);
		solver.A.initApply();
	}
};

template <> struct AX_MKLPDSSSystem<AX_AcousticRealLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_AcousticRealLinear *analysis)
	: AX_MKLPDSSSystemData(analysis->configuration.mklpdss)
	{
		assembler.A.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
		_fillAssembler(this, analysis->configuration.acoustic_pressure, 1);

		solver.A.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
		_fillSolver(this, analysis->configuration.acoustic_pressure, 2);
	}
};

template <> struct AX_MKLPDSSSystem<AX_AcousticComplexLinear>: public AX_MKLPDSSSystemData<double, std::complex<double> > {

	AX_MKLPDSSSystem(AX_AcousticComplexLinear *analysis)
	: AX_MKLPDSSSystemData(analysis->configuration.mklpdss)
	{
		assembler.A.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
		_fillAssembler(this, analysis->configuration.acoustic_pressure, 1);

		if (analysis->configuration.impedance.size()) {
			solver.A.type = Matrix_Type::COMPLEX_SYMMETRIC;
		} else {
			solver.A.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
		}

		_fillSolver(this, analysis->configuration.acoustic_pressure, 1);
	}
};

void setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet, const DOFsDistribution &distribution);
void setDirichlet(Matrix_Distributed<Matrix_CSR, std::complex<double>> &A, Vector_Distributed<Vector_Dense, std::complex<double>> &b, const Vector_Sparse<std::complex<double>> &dirichlet, const DOFsDistribution &distribution);

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_ */
