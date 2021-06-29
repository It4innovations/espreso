
#ifndef SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_STEADYSTATE_LINEAR_H_

#include "analysis/scheme/steadystate.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsLoadStepConfiguration;

class AX_StructuralMechanicsSteadyStateLinear {

public:
	AX_StructuralMechanicsSteadyStateLinear(StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);

	StructuralMechanicsGlobalSettings &gsettings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	Matrix_Base<double> *K;
	Vector_Base<double> *f;

	AX_SteadyState scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_STEADYSTATE_LINEAR_H_ */
