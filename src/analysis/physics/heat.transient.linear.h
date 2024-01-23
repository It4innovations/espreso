
#ifndef SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_
#define SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/heattransfer.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatTransientLinear: public Physics {

public:
	HeatTransientLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
	~HeatTransientLinear();

	void analyze();
	void run(step::Step &step);

	step::Time time;
	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	HeatTransfer assembler;

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *f, *x, *dirichlet;
	Vector_Base<double> *U, *dU, *V, *X, *Y, *dTK, *dTM;

	SparseMatrixBuilder<double> *builder;
	LinearSystemSolver<double> *solver;

protected:
	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}



#endif /* SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_ */
