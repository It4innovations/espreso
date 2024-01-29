
#ifndef SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsSteadyStateLinear: public Physics {

public:
	StructuralMechanicsSteadyStateLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanicsSteadyStateLinear();

	void analyze(step::Step &step);
	void run(step::Step &step);

	step::Time time;
	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	StructuralMechanics assembler;

	Matrix_Base<double> *K;
	Vector_Base<double> *f, *x, *dirichlet;

	SparseMatrixBuilder<double> *builder;
	LinearSystemSolver<double> *solver;

protected:
	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_ */
