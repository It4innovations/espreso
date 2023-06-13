
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct NonLinearSolverConfiguration;

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

	Matrix_Base<double> *K;
	Vector_Base<double> *U, *R, *f, *x, *dirichlet;

	LinearSystem<double> *system;

protected:
	bool checkTemp(step::Step &step);

	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_ */
