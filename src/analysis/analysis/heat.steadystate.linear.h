
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class AX_HeatSteadyStateLinear: public Analysis {

public:
	AX_HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~AX_HeatSteadyStateLinear();

	void analyze();
	void run(step::Step &step);

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	AX_HeatTransfer assembler;
	AX_SteadyState scheme;

	AX_LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_ */
