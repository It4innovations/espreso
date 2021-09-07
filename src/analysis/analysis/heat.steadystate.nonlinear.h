
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/solver/newtonraphson.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferGlobalSettings;
struct HeatTransferLoadStepConfiguration;

class AX_HeatSteadyStateNonLinear: public Analysis {

public:
	AX_HeatSteadyStateNonLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);
	~AX_HeatSteadyStateNonLinear();

	void init();
	void run(step::Step &step);

	HeatTransferGlobalSettings &gsettings;
	HeatTransferLoadStepConfiguration &configuration;

	AX_HeatTransfer assembler;
	AX_NewtonRaphson solver;
	AX_SteadyState scheme;

	AX_LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_ */
