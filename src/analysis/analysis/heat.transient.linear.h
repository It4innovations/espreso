
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_LINEAR_H_

#include "analysis/scheme/transient.firstorderimplicit.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct HeatTransferGlobalSettings;
struct HeatTransferLoadStepConfiguration;

class AX_HeatTransientLinear {

public:
	AX_HeatTransientLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);

	HeatTransferGlobalSettings &gsettings;
	HeatTransferLoadStepConfiguration &configuration;

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *f;

	AX_TransientFirstOrderImplicit scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_TRANSIENT_LINEAR_H_ */
