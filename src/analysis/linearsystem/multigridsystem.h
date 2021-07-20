
#ifndef SRC_ANALYSIS_LINEARSOLVER_MULTIGRIDSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_MULTIGRIDSOLVER_H_

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

#include "config/ecf/linearsolver/hypre/hypre.h"

namespace espreso {

template <typename T>
struct AX_MultigridSystem: AX_LinearSystem<T> {

	AX_MultigridSystem(HYPREConfiguration &configuration) {}

	void setMapping(Matrix_Base<T> *A) const
	{

	}

	void setMapping(Vector_Base<T> *x) const
	{

	}

	void init(AX_AcousticRealLinear *analysis)
	{

	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{

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

	void update(AX_Acoustic &assembler)
	{

	}

	void update(AX_HeatTransfer &assembler)
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

	bool solve()
	{
		printf("Multigrid solve\n");
		return false;
	}

//	Matrix_FETI<Matrix_CSR, T> K;
//
//	Matrix_FETI<Matrix_Dense, T> N1, N2;
//	Matrix_FETI<Matrix_IJV, T> B1, B0;
//
//	Vector_FETI<Vector_Dense, T> B1c, B1duplication, B1gap;

//	FETISystemSolver solver;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_MULTIGRIDSOLVER_H_ */
