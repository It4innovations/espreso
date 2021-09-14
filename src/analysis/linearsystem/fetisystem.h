
#ifndef SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_

#include "linearsystem.h"
#include "analysis/analysis/heat.steadystate.linear.h"
#include "analysis/analysis/heat.steadystate.nonlinear.h"
#include "analysis/analysis/acoustic.real.linear.h"
//#include "analysis/analysis/heat.transient.linear.h"
//#include "analysis/analysis/heat.transient.nonlinear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.complex.linear.h"
//#include "analysis/analysis/structuralmechanics.harmonic.real.linear.h"
//#include "analysis/analysis/structuralmechanics.steadystate.linear.h"
//#include "analysis/analysis/structuralmechanics.transient.linear.h"

#include "analysis/composer/nodes.uniform.feti.h"
#include "config/ecf/linearsolver/feti.h"

#include "math2/generalization/matrix_feti.h"

namespace espreso {

template <typename T>
struct AX_FETISystem: AX_LinearSystem<T> {

	AX_FETISystem(FETIConfiguration &configuration): configuration(configuration) {}

	void setMapping(Matrix_Base<T> *A) const
	{
		pattern.setMap(dynamic_cast<Matrix_FETI<Matrix_CSR, T>*>(A));
	}

	void setMapping(Vector_Base<T> *x) const
	{
		pattern.setMap(dynamic_cast<Vector_FETI<Vector_Dense, T>*>(x));
	}

	void init(AX_AcousticRealLinear *analysis)
	{

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

//		analysis->assembler.initDirichlet(dirichlet);
//		this->assembler.dirichlet = this->solver.dirichlet = &dirichlet;

		initKernels(*this, analysis->assembler);
	}

	void init(AX_HeatSteadyStateNonLinear *analysis)
	{
		pattern.set(1, analysis->assembler.matrixType());
		pattern.fill(A);
		pattern.fill(b);
		pattern.fill(x);
		this->assembler.A = this->solver.A = &A;
		this->assembler.b = this->solver.b = &b;
		this->assembler.x = this->solver.x = &x;

//		analysis->assembler.initDirichlet(dirichlet);
//		this->assembler.dirichlet = this->solver.dirichlet = &dirichlet;

		initKernels(*this, analysis->assembler);
	}

//	void init(AX_HeatSteadyStateNonLinear *analysis)
//	{
//		pattern.set(1, analysis->assembler.matrixType());
//	}
//
//	void init(AX_HeatTransientLinear *analysis)
//	{
//		pattern.set(1, analysis->assembler.matrixType());
//	}
//
//	void init(AX_HeatTransientNonLinear *analysis)
//	{
//		pattern.set(1, analysis->assembler.matrixType());
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
//		if (solver.A.touched || solver.b.touched || solver.dirichlet.touched) {
			updateKernels(*this, assembler);
//		}
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
//
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
		printf("FETI solve\n");
		return false;
	}

	FETIConfiguration &configuration;
	UniformNodesFETIPattern pattern;

	Matrix_FETI<Matrix_CSR, T> A;
	Vector_FETI<Vector_Dense, T> x, b;
	Vector_Sparse<T> dirichlet;

	Matrix_FETI<Matrix_CSR, T> RegMat;
	Matrix_FETI<Matrix_Dense, T> N1, N2;
	Matrix_FETI<Matrix_IJV, T> B1, B0;

	Vector_FETI<Vector_Dense, T> B1c, B1duplication;

//	FETISystemSolver solver;
};

void initGluing(AX_FETISystem<double> &solver);

void initKernels(AX_FETISystem<double> &solver, AX_HeatTransfer &assembler);
void updateKernels(AX_FETISystem<double> &solver, AX_HeatTransfer &assembler);

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_FETISOLVER_H_ */
