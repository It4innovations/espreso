
#ifndef SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_

#include "math2/math2.h"

namespace espreso {

class AX_HeatSteadyStateLinear;
class AX_HeatSteadyStateNonLinear;
class AX_HeatTransientLinear;
class AX_HeatTransientNonLinear;
class AX_StructuralMechanicsHarmonicComplexLinear;
class AX_StructuralMechanicsHarmonicRealLinear;
class AX_StructuralMechanicsSteadyStateLinear;
class AX_StructuralMechanicsTransientLinear;

template <typename T> class DirectSolver;
template <typename T> class FETISolver;
template <typename T> class MultigridSolver;

struct LinearSolver {

	virtual ~LinearSolver() {}

	virtual void init(AX_HeatSteadyStateLinear    *analysis) =0;
	virtual void init(AX_HeatSteadyStateNonLinear *analysis) =0;
	virtual void init(AX_HeatTransientLinear      *analysis) =0;
	virtual void init(AX_HeatTransientNonLinear   *analysis) =0;
	virtual void init(AX_StructuralMechanicsHarmonicComplexLinear *analysis) =0;
	virtual void init(AX_StructuralMechanicsHarmonicRealLinear    *analysis) =0;
	virtual void init(AX_StructuralMechanicsSteadyStateLinear     *analysis) =0;
	virtual void init(AX_StructuralMechanicsTransientLinear       *analysis) =0;

	virtual void prepare(AX_HeatSteadyStateLinear    *analysis) =0;
	virtual void prepare(AX_HeatSteadyStateNonLinear *analysis) =0;
	virtual void prepare(AX_HeatTransientLinear      *analysis) =0;
	virtual void prepare(AX_HeatTransientNonLinear   *analysis) =0;
	virtual void prepare(AX_StructuralMechanicsHarmonicComplexLinear *analysis) =0;
	virtual void prepare(AX_StructuralMechanicsHarmonicRealLinear    *analysis) =0;
	virtual void prepare(AX_StructuralMechanicsSteadyStateLinear     *analysis) =0;
	virtual void prepare(AX_StructuralMechanicsTransientLinear       *analysis) =0;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_ */
