
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_

#include "analysis/assembler/assembler.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/mode/linear.h"
#include "analysis/linearsolver/linearsolver.h"

namespace espreso {

struct HeatTransferGlobalSettings;
struct HeatTransferLoadStepConfiguration;

class AX_HeatSteadyStateLinear {

public:
	AX_HeatSteadyStateLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);

	void solve();

	HeatTransferGlobalSettings &gsettings;
	HeatTransferLoadStepConfiguration &configuration;

	Assembler assembler;
	AX_SteadyState scheme;
	AX_Linear mode;
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_ */
