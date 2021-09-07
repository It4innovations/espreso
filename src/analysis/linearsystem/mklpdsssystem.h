
#ifndef SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/acoustic.real.linear.h"
//#include "analysis/analysis/heat.steadystate.nonlinear.h"
//#include "analysis/analysis/heat.transient.linear.h"
//#include "analysis/analysis/heat.transient.nonlinear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.complex.linear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.real.linear.h"
//#include "analysis/analysis/structuralmechanics.steadystate.linear.h"
//#include "analysis/analysis/structuralmechanics.transient.linear.h"

#include "analysis/composer/nodes.uniform.distributed.h"
#include "basis/utilities/sysutils.h"
#include "config/ecf/linearsolver/mklpdss.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "math2/generalization/matrix_distributed.h"
#include "wrappers/mklpdss/w.mkl.pdss.h"

namespace espreso {

template <typename T>
struct AX_MKLPDSSSystem: public AX_LinearSystem<T> {

	AX_MKLPDSSSystem(MKLPDSSConfiguration &configuration): mklpdss(configuration), set(true) {}

	void setMapping(Matrix_Base<T> *A) const
	{
		assembler.pattern.setMap(dynamic_cast<Matrix_Distributed<Matrix_CSR, T>*>(A));
	}

	void setMapping(Vector_Base<T> *x) const
	{
		assembler.pattern.setMap(dynamic_cast<Vector_Distributed<Vector_Dense, T>*>(x));
	}

	void info() const
	{
		mklpdss.info(solver.A);
	}

	void init(AX_AcousticRealLinear *analysis)
	{
		assembler.A.type = analysis->assembler.matrixType();
		assembler.pattern.set(1);
		assembler.pattern.fill(assembler.A);
		assembler.pattern.fill(assembler.b);
		assembler.pattern.fill(assembler.x);
		AX_LinearSystem<T>::assembler.A = &assembler.A;
		AX_LinearSystem<T>::assembler.b = &assembler.b;
		AX_LinearSystem<T>::assembler.x = &assembler.x;
		AX_LinearSystem<T>::assembler.dirichlet = &assembler.dirichlet;

		solver.A.type = analysis->assembler.matrixType();
		solver.pattern.set(2);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		AX_LinearSystem<T>::solver.A = &solver.A;
		AX_LinearSystem<T>::solver.b = &solver.b;
		AX_LinearSystem<T>::solver.x = &solver.x;
		AX_LinearSystem<T>::solver.dirichlet = &solver.dirichlet;

		analysis->assembler.initDirichlet(assembler.dirichlet);
		math::multiplyPattern(solver.dirichlet, assembler.dirichlet, 1, 2);
	}

	template<typename HeatSteadyState>
	void _initHeat(HeatSteadyState *analysis)
	{
		assembler.A.type = solver.A.type = analysis->assembler.matrixType();
		assembler.pattern.set(1);
		assembler.pattern.fill(solver.A);
		assembler.pattern.fill(solver.b);
		assembler.pattern.fill(solver.x);
		AX_LinearSystem<T>::assembler.A = AX_LinearSystem<T>::solver.A = &solver.A;
		AX_LinearSystem<T>::assembler.b = AX_LinearSystem<T>::solver.b = &solver.b;
		AX_LinearSystem<T>::assembler.x = AX_LinearSystem<T>::solver.x = &solver.x;

		analysis->assembler.initDirichlet(solver.dirichlet);
		AX_LinearSystem<T>::assembler.dirichlet = AX_LinearSystem<T>::solver.dirichlet = &solver.dirichlet;
	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{
		_initHeat(analysis);
	}

	void init(AX_HeatSteadyStateNonLinear *analysis)
	{
		_initHeat(analysis);
	}

//	void init(AX_HeatSteadyStateNonLinear *analysis)
//	{
//
//	}
//
//	void init(AX_HeatTransientLinear *analysis)
//	{
//
//	}
//
//	void init(AX_HeatTransientNonLinear *analysis)
//	{
//
//	}
//
//	void init(AX_StructuralMechanicsHarmonicComplexLinear *analysis)
//	{
//
//	}
//
//	void init(AX_StructuralMechanicsHarmonicRealLinear *analysis)
//	{
//
//	}
//
//	void init(AX_StructuralMechanicsSteadyStateLinear *analysis)
//	{
//
//	}
//
//	void init(AX_StructuralMechanicsTransientLinear *analysis)
//	{
//
//	}

	void _update(step::Step &step)
	{
		if (set) {
			mklpdss.set(solver.A);
			set = false;
		}
		if (solver.A.touched || solver.b.touched || solver.dirichlet.touched) {
			setDirichlet(solver.A, solver.b, solver.dirichlet);
			mklpdss.update(solver.A);
		}

		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: system/{A, b, dirichlet}\n");
			math::store(solver.A, utils::filename(utils::debugDirectory(step) + "/system", "A").c_str());
			math::store(solver.b, utils::filename(utils::debugDirectory(step) + "/system", "b").c_str());
			math::store(solver.dirichlet, utils::filename(utils::debugDirectory(step) + "/system", "dirichlet").c_str());
		}
	}

	void update(step::Step &step, AX_Acoustic &assembler)
	{
		_update(step);
	}

	void update(step::Step &step, AX_HeatTransfer &assembler)
	{
		_update(step);
	}

//	void prepare(AX_HeatSteadyStateNonLinear *analysis)
//	{
//
//	}
//
//	void prepare(AX_HeatTransientLinear *analysis)
//	{
//
//	}
//
//	void prepare(AX_HeatTransientNonLinear *analysis)
//	{
//
//	}
//
//	void prepare(AX_StructuralMechanicsHarmonicComplexLinear *analysis)
//	{
//
//	}
//	void prepare(AX_StructuralMechanicsHarmonicRealLinear *analysis)
//	{
//
//	}
//
//	void prepare(AX_StructuralMechanicsSteadyStateLinear *analysis)
//	{
//
//	}
//
//	void prepare(AX_StructuralMechanicsTransientLinear *analysis)
//	{
//
//	}

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

	struct Data {
		UniformNodesDistributedPattern pattern;

		Matrix_Distributed<Matrix_CSR, T> A;
		Vector_Distributed<Vector_Dense, T> x, b;
		Vector_Sparse<T> dirichlet;
	} assembler, solver;

	MKLPDSS<T> mklpdss;
	bool set;
};

void setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet);

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_ */
