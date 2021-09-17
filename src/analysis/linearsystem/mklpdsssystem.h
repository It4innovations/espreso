
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
#include "math2/utils/dofs_distribution.h"
#include "math2/utils/utils_distributed.h"
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
			setDirichlet(solver.A, solver.b, solver.dirichlet.cluster, solver.A.distribution);
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
	};
	
	Data<Assembler> assembler;
	Data<Solver> solver;

	MKLPDSS<Solver> mklpdss;
};

template <typename Analysis> struct AX_MKLPDSSSystem {};

inline void _initDirect(AX_MKLPDSSSystemData<double, double> &system, std::map<std::string, ECFExpression> &dirichlet)
{
	system.assembler.pattern.set(dirichlet, 1, system.solver.A.distribution, system.solver.A.synchronization);
	system.assembler.pattern.fill(system.solver.A);
	system.assembler.pattern.fill(system.solver.b);
	system.assembler.pattern.fill(system.solver.x);
	system.assembler.pattern.fill(system.solver.dirichlet);

	system.AX_LinearSystem::assembler.A = system.AX_LinearSystem::solver.A = &system.solver.A;
	system.AX_LinearSystem::assembler.b = system.AX_LinearSystem::solver.b = &system.solver.b;
	system.AX_LinearSystem::assembler.x = system.AX_LinearSystem::solver.x = &system.solver.x;
	system.AX_LinearSystem::assembler.dirichlet = system.AX_LinearSystem::solver.dirichlet = &system.solver.dirichlet;
}

template <> struct AX_MKLPDSSSystem<AX_HeatSteadyStateLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_HeatSteadyStateLinear *analysis, MKLPDSSConfiguration &configuration)
	: AX_MKLPDSSSystemData(configuration)
	{
		solver.A.type = analysis->assembler.matrixType();
		_initDirect(*this, analysis->configuration.temperature);
	}
};

template <> struct AX_MKLPDSSSystem<AX_HeatSteadyStateNonLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_HeatSteadyStateNonLinear *analysis, MKLPDSSConfiguration &configuration)
	: AX_MKLPDSSSystemData(configuration)
	{
		solver.A.type = analysis->assembler.matrixType();
		_initDirect(*this, analysis->configuration.temperature);
	}
};

template <> struct AX_MKLPDSSSystem<AX_AcousticRealLinear>: public AX_MKLPDSSSystemData<double, double> {

	AX_MKLPDSSSystem(AX_AcousticRealLinear *analysis, MKLPDSSConfiguration &configuration)
	: AX_MKLPDSSSystemData(configuration)
	{
		assembler.A.type = analysis->assembler.matrixType();
		assembler.pattern.set(analysis->configuration.acoustic_pressure, 1);
		assembler.pattern.fill(assembler.A);
		assembler.pattern.fill(assembler.b);
		assembler.pattern.fill(assembler.x);
		assembler.pattern.fill(assembler.dirichlet);
		this->AX_LinearSystem::assembler.A = &assembler.A;
		this->AX_LinearSystem::assembler.b = &assembler.b;
		this->AX_LinearSystem::assembler.x = &assembler.x;
		this->AX_LinearSystem::assembler.dirichlet = &assembler.dirichlet;

		solver.A.type = analysis->assembler.matrixType();
		solver.pattern.set(analysis->configuration.acoustic_pressure, 2, solver.A.distribution, solver.A.synchronization);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		solver.pattern.fill(solver.dirichlet);
		this->AX_LinearSystem::solver.A = &solver.A;
		this->AX_LinearSystem::solver.b = &solver.b;
		this->AX_LinearSystem::solver.x = &solver.x;
		this->AX_LinearSystem::solver.dirichlet = &solver.dirichlet;
	}
};

template <> struct AX_MKLPDSSSystem<AX_AcousticComplexLinear>: public AX_MKLPDSSSystemData<double, std::complex<double> > {

	AX_MKLPDSSSystem(AX_AcousticComplexLinear *analysis, MKLPDSSConfiguration &configuration)
	: AX_MKLPDSSSystemData(configuration)
	{
		assembler.A.type = analysis->assembler.matrixType();
		assembler.pattern.set(analysis->configuration.acoustic_pressure, 1);
		assembler.pattern.fill(assembler.A);
		assembler.pattern.fill(assembler.b);
		assembler.pattern.fill(assembler.x);
		assembler.pattern.fill(assembler.dirichlet);
		this->AX_LinearSystem::assembler.A = &assembler.A;
		this->AX_LinearSystem::assembler.b = &assembler.b;
		this->AX_LinearSystem::assembler.x = &assembler.x;
		this->AX_LinearSystem::assembler.dirichlet = &assembler.dirichlet;

		solver.A.type = analysis->assembler.matrixType();
		solver.pattern.set(analysis->configuration.acoustic_pressure, 1, solver.A.distribution, solver.A.synchronization);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		solver.pattern.fill(solver.dirichlet);
		this->AX_LinearSystem::solver.A = &solver.A;
		this->AX_LinearSystem::solver.b = &solver.b;
		this->AX_LinearSystem::solver.x = &solver.x;
		this->AX_LinearSystem::solver.dirichlet = &solver.dirichlet;
	}
};

void setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet, const DOFsDistribution &distribution);
void setDirichlet(Matrix_Distributed<Matrix_CSR, std::complex<double>> &A, Vector_Distributed<Vector_Dense, std::complex<double>> &b, const Vector_Sparse<std::complex<double>> &dirichlet, const DOFsDistribution &distribution);

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_ */
