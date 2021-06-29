
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_NONLINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_NONLINEAR_H_

#include "analysis/scheme/transient.firstorderimplicit.h"
#include "analysis/mode/newtonraphson.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct HeatTransferGlobalSettings;
struct HeatTransferLoadStepConfiguration;

class AX_HeatTransientNonLinear {

public:
	AX_HeatTransientNonLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);

	HeatTransferGlobalSettings &gsettings;
	HeatTransferLoadStepConfiguration &configuration;

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *f;

	AX_TransientFirstOrderImplicit scheme;
	AX_NewtonRaphson mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_NONLINEAR_H_ */
