
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"
#include "analysis/nonlinearity/newtonraphson.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatSteadyStateNonLinear: public Analysis {

public:
	HeatSteadyStateNonLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~HeatSteadyStateNonLinear();

	void analyze();
	void run(step::Step &step);

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	HeatTransfer assembler;
	NewtonRaphson solver;
	SteadyState scheme;

	LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_ */
