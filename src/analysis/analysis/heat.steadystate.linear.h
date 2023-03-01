
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatSteadyStateLinear: public Analysis {

public:
	HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~HeatSteadyStateLinear();

	void analyze();
	void run(step::Step &step);

	void dryrun();

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	HeatTransfer assembler;
	SteadyState scheme;

	LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_ */
