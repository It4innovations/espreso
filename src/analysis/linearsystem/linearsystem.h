
#ifndef SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_

#include "analysis/composer/elementmapping.h"
#include "math2/math2.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

namespace step { struct Step; }

class AX_HeatSteadyStateLinear;
class AX_HeatSteadyStateNonLinear;
class AX_AcousticRealLinear;
//class AX_HeatTransientLinear;
//class AX_HeatTransientNonLinear;
//class AX_StructuralMechanicsHarmonicComplexLinear;
//class AX_StructuralMechanicsHarmonicRealLinear;
//class AX_StructuralMechanicsSteadyStateLinear;
//class AX_StructuralMechanicsTransientLinear;

class AX_HeatTransfer;
class AX_Acoustic;

template <typename T>
struct AX_LinearSystem {

	struct Template {
		Matrix_Base<T> *A;
		Vector_Base<T> *x, *b;
		Vector_Sparse<T> *dirichlet;
	};

	virtual void setMapping(Matrix_Base<T> *A) const =0;
	virtual void setMapping(Vector_Base<T> *x) const =0;

	virtual ~AX_LinearSystem() {}

	virtual void info() const {};

	virtual void init(AX_AcousticRealLinear       *analysis) =0;
	virtual void init(AX_HeatSteadyStateLinear    *analysis) =0;
	virtual void init(AX_HeatSteadyStateNonLinear *analysis) =0;
//	virtual void init(AX_HeatSteadyStateNonLinear *analysis) =0;
//	virtual void init(AX_HeatTransientLinear      *analysis) =0;
//	virtual void init(AX_HeatTransientNonLinear   *analysis) =0;
//	virtual void init(AX_StructuralMechanicsHarmonicComplexLinear *analysis) =0;
//	virtual void init(AX_StructuralMechanicsHarmonicRealLinear    *analysis) =0;
//	virtual void init(AX_StructuralMechanicsSteadyStateLinear     *analysis) =0;
//	virtual void init(AX_StructuralMechanicsTransientLinear       *analysis) =0;

	virtual void update(step::Step &step, AX_Acoustic     &assembler) =0;
	virtual void update(step::Step &step, AX_HeatTransfer &assembler) =0;
//	virtual void update(AX_HeatSteadyStateNonLinear *analysis, bool A, bool b) =0;
//	virtual void update(AX_HeatTransientLinear      *analysis, bool A, bool b) =0;
//	virtual void update(AX_HeatTransientNonLinear   *analysis, bool A, bool b) =0;
//	virtual void update(AX_StructuralMechanicsHarmonicComplexLinear *analysis, bool A, bool b) =0;
//	virtual void update(AX_StructuralMechanicsHarmonicRealLinear    *analysis, bool A, bool b) =0;
//	virtual void update(AX_StructuralMechanicsSteadyStateLinear     *analysis, bool A, bool b) =0;
//	virtual void update(AX_StructuralMechanicsTransientLinear       *analysis, bool A, bool b) =0;

	virtual bool solve(step::Step &step) =0;

	Template assembler, solver;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_ */
