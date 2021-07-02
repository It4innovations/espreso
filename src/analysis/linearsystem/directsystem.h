
#ifndef SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
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

	ElementMapping<T> mapping(const Matrix_Base<T> *A) const
	{
		return ElementMapping<T>();
	}

	ElementMapping<T> mapping(const Vector_Base<T> *x) const
	{
		return ElementMapping<T>();
	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{
		pattern.set(1, analysis->assembler.matrixType());
		pattern.fill(A);
		pattern.fill(b);
		pattern.fill(x);
		this->assembler.A = this->solver.A = &A;
		this->assembler.b = this->solver.b = &b;
		this->assembler.x = this->solver.x = &x;
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

	void update(AX_HeatTransfer &assembler, bool A, bool b, bool dirichlet)
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
		printf("Direct solve\n");
		return false;
	}

	UniformNodesDistributedPattern pattern;

	Matrix_Distributed<Matrix_CSR, T> A;
	Vector_Distributed<Vector_Dense, T> x, b;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_DIRECTSOLVER_H_ */
