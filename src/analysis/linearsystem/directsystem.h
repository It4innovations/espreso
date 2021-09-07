
#ifndef SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/acoustic.real.linear.h"
//#include "analysis/analysis/heat.steadystate.nonlinear.h"
//#include "analysis/analysis/heat.transient.linear.h"
//#include "analysis/analysis/heat.transient.nonlinear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.complex.linear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.real.linear.h"
//#include "analysis/analysis/structuralmechanics.steadystate.linear.h"
//#include "analysis/analysis/structuralmechanics.transient.linear.h"

#include "analysis/composer/nodes.uniform.distributed.h"

#include "config/ecf/linearsolver/mklpdss.h"
#include "config/ecf/linearsolver/pardiso.h"
#include "config/ecf/linearsolver/superlu.h"
#include "config/ecf/linearsolver/wsmp.h"

#include "math2/generalization/matrix_distributed.h"

namespace espreso {

template <typename T>
struct AX_DirectSystem: AX_LinearSystem<T> {

	AX_DirectSystem(MKLPDSSConfiguration &configuration) {}
	AX_DirectSystem(PARDISOConfiguration &configuration) {}
	AX_DirectSystem(SuperLUConfiguration &configuration) {}
	AX_DirectSystem(WSMPConfiguration &configuration) {}

	void setMapping(Matrix_Base<T> *A) const
	{

	}

	void setMapping(Vector_Base<T> *x) const
	{

	}

	void init(AX_AcousticRealLinear *analysis)
	{
		assembler.pattern.set(1);
		assembler.pattern.fill(assembler.A);
		assembler.pattern.fill(assembler.b);
		assembler.pattern.fill(assembler.x);
		AX_LinearSystem<T>::assembler.A = &assembler.A;
		AX_LinearSystem<T>::assembler.b = &assembler.b;
		AX_LinearSystem<T>::assembler.x = &assembler.x;

		solver.pattern.set(2);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		AX_LinearSystem<T>::solver.A = &solver.A;
		AX_LinearSystem<T>::solver.b = &solver.b;
		AX_LinearSystem<T>::solver.x = &solver.x;
	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{
		solver.pattern.set(1);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		AX_LinearSystem<T>::assembler.A = AX_LinearSystem<T>::solver.A = &solver.A;
		AX_LinearSystem<T>::assembler.b = AX_LinearSystem<T>::solver.b = &solver.b;
		AX_LinearSystem<T>::assembler.x = AX_LinearSystem<T>::solver.x = &solver.x;
	}

	void init(AX_HeatSteadyStateNonLinear *analysis)
	{
		solver.pattern.set(1);
		solver.pattern.fill(solver.A);
		solver.pattern.fill(solver.b);
		solver.pattern.fill(solver.x);
		AX_LinearSystem<T>::assembler.A = AX_LinearSystem<T>::solver.A = &solver.A;
		AX_LinearSystem<T>::assembler.b = AX_LinearSystem<T>::solver.b = &solver.b;
		AX_LinearSystem<T>::assembler.x = AX_LinearSystem<T>::solver.x = &solver.x;
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

	void update(step::Step &step, AX_Acoustic &assembler)
	{

	}

	void update(step::Step &step, AX_HeatTransfer &assembler)
	{

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
		printf("Direct solve\n");
		return false;
	}

	struct Data {
		UniformNodesDistributedPattern pattern;

		Matrix_Distributed<Matrix_CSR, T> A;
		Vector_Distributed<Vector_Dense, T> x, b;
	} assembler, solver;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_ */
