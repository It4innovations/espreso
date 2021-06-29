
#ifndef SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_COMPLEX_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_COMPLEX_LINEAR_H_

#include "analysis/scheme/harmonic.complex.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

#include <complex>

namespace espreso {

struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsLoadStepConfiguration;

class AX_StructuralMechanicsHarmonicComplexLinear {

public:
	AX_StructuralMechanicsHarmonicComplexLinear(StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);

	StructuralMechanicsGlobalSettings &gsettings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	Matrix_Base<std::complex<double> > *K;
	Vector_Base<std::complex<double> > *f;

	AX_HarmonicComplex scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_STRUCTURALMECHANICS_HARMONIC_COMPLEX_LINEAR_H_ */
