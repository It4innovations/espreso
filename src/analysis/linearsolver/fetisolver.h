
#ifndef SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_

#include "linearsolver.h"

#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/heat.transient.linear.h"
#include "analysis/analysis/heat.transient.nonlinear.h"
#include "analysis/analysis/structuralmechanics.harmonic.complex.linear.h"
#include "analysis/analysis/structuralmechanics.harmonic.real.linear.h"
#include "analysis/analysis/structuralmechanics.steadystate.linear.h"
#include "analysis/analysis/structuralmechanics.transient.linear.h"

#include "analysis/composer/nodes.uniform.feti.h"
#include "config/ecf/linearsolver/feti.h"

#include "math2/generalization/matrix_feti.h"

namespace espreso {

template <typename T>
struct FETISolver: LinearSolver {

	FETISolver(FETIConfiguration &configuration) {}

	Matrix_FETI<Matrix_CSR, T>* getMatrixPattern()
	{
		return NULL;
	}

	Vector_FETI<Vector_Dense, T>* getVectorPattern()
	{
		return NULL;
	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{
		printf("init FETI pattern\n");
		switch (analysis->assembler.matrixType()) {
		case Matrix_Type::REAL_SYMMETRIC_INDEFINITE: pattern.initUpper(); break;
		case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE: pattern.initUpper(); break;
		case Matrix_Type::REAL_UNSYMMETRIC: pattern.initFull(); break;
		}
	}

	void init(AX_HeatSteadyStateNonLinear *analysis)
	{

	}

	void init(AX_HeatTransientLinear *analysis)
	{

	}

	void init(AX_HeatTransientNonLinear *analysis)
	{

	}

	void init(AX_StructuralMechanicsHarmonicComplexLinear *analysis)
	{

	}

	void init(AX_StructuralMechanicsHarmonicRealLinear *analysis)
	{

	}

	void init(AX_StructuralMechanicsSteadyStateLinear *analysis)
	{

	}

	void init(AX_StructuralMechanicsTransientLinear *analysis)
	{

	}

	void prepare(AX_HeatSteadyStateLinear *analysis)
	{

	}

	void prepare(AX_HeatSteadyStateNonLinear *analysis)
	{

	}

	void prepare(AX_HeatTransientLinear *analysis)
	{

	}

	void prepare(AX_HeatTransientNonLinear *analysis)
	{

	}

	void prepare(AX_StructuralMechanicsHarmonicComplexLinear *analysis)
	{

	}

	void prepare(AX_StructuralMechanicsHarmonicRealLinear *analysis)
	{

	}

	void prepare(AX_StructuralMechanicsSteadyStateLinear *analysis)
	{

	}

	void prepare(AX_StructuralMechanicsTransientLinear *analysis)
	{

	}

	UniformNodesFETIPattern pattern;

	Matrix_FETI<Matrix_CSR, T> K;

	Matrix_FETI<Matrix_Dense, T> N1, N2;
	Matrix_FETI<Matrix_IJV, T> B1, B0;

	Vector_FETI<Vector_Dense, T> B1c, B1duplication, B1gap;

//	FETISystemSolver solver;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_ */
