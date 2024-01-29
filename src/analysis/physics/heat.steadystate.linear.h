
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/heattransfer.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatSteadyStateLinear: public Physics {

public:
	HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~HeatSteadyStateLinear();

	void analyze(step::Step &step);
	void run(step::Step &step);

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	HeatTransfer assembler;

	Matrix_Base<double> *K;
	Vector_Base<double> *f, *x, *dirichlet;

	SparseMatrixBuilder<double> *builder;
	LinearSystemSolver<double> *solver;

protected:
	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_ */
