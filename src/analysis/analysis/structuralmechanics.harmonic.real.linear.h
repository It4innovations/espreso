
#ifndef SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_

#include "analysis/scheme/harmonic.real.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsLoadStepConfiguration;

class AX_StructuralMechanicsHarmonicRealLinear {

public:
	AX_StructuralMechanicsHarmonicRealLinear(StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);

	StructuralMechanicsGlobalSettings &gsettings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	Matrix_Base<double> *K;
	Vector_Base<double> *f;

	AX_HarmonicReal scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_ */
