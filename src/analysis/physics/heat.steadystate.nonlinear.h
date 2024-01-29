
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/heattransfer.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct NonLinearSolverConfiguration;

class HeatSteadyStateNonLinear: public Physics {

public:
	HeatSteadyStateNonLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~HeatSteadyStateNonLinear();

	void analyze(step::Step &step);
	void run(step::Step &step);

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	HeatTransfer assembler;

	Matrix_Base<double> *K;
	Vector_Base<double> *U, *R, *f, *x, *dirichlet;

	SparseMatrixBuilder<double> *builder;
	LinearSystemSolver<double> *solver;

protected:
	bool checkTemp(step::Step &step);

	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_NONLINEAR_H_ */
