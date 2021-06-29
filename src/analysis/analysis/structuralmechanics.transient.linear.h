
#ifndef SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_

#include "analysis/scheme/transient.secondorderimplicit.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsLoadStepConfiguration;

class AX_StructuralMechanicsTransientLinear {

public:
	AX_StructuralMechanicsTransientLinear(StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);

	StructuralMechanicsGlobalSettings &gsettings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *f;

	AX_TransientSecondOrderImplicit scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_ */
