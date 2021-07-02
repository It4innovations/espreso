
#ifndef SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_

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
#include "basis/utilities/sysutils.h"
#include "config/ecf/linearsolver/mklpdss.h"
#include "esinfo/ecfinfo.h"
#include "math2/generalization/matrix_distributed.h"
#include "wrappers/mklpdss/w.mkl.pdss.h"

namespace espreso {

template <typename T>
struct AX_MKLPDSSSystem: public AX_LinearSystem<T> {

	AX_MKLPDSSSystem(MKLPDSSConfiguration &configuration): mklpdss(configuration) {}

	ElementMapping<T> mapping(const Matrix_Base<T> *A) const
	{
		return pattern.map(dynamic_cast<const Matrix_Distributed<Matrix_CSR, T>*>(A));
	}

	ElementMapping<T> mapping(const Vector_Base<T> *x) const
	{
		return pattern.map(dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(x));
	}

	void init(AX_HeatSteadyStateLinear *analysis)
	{
		pattern.set(1, A.cluster.type = analysis->assembler.matrixType());
		pattern.fill(A);
		pattern.fill(b);
		pattern.fill(x);
		this->assembler.A = this->solver.A = &A;
		this->assembler.b = this->solver.b = &b;
		this->assembler.x = this->solver.x = &x;
		analysis->assembler.fillDirichletIndices(this->dirichlet);
		mklpdss.set(A);
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
		if (A || b || dirichlet) {
			setDirichlet(this->A, this->b, this->dirichlet);
			mklpdss.update(this->A);
		}

		if (info::ecf->output.print_matrices) {
			math::store(this->A, utils::filename(utils::debugDirectory() + "/system", "A").c_str());
			math::store(this->b, utils::filename(utils::debugDirectory() + "/system", "b").c_str());
		}
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
		return mklpdss.solve(b, x);
	}

	UniformNodesDistributedPattern pattern;

	Matrix_Distributed<Matrix_CSR, T> A;
	Vector_Distributed<Vector_Dense, T> x, b;

	MKLPDSS<T> mklpdss;
};

void setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet);

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_MKLPDSSSOLVER_H_ */
