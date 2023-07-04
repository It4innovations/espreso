
#ifndef SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_
#define SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/module/structuralmechanics.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsSteadyStateNonLinear: public Physics {

public:
	StructuralMechanicsSteadyStateNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanicsSteadyStateNonLinear();

	void analyze();
	void run(step::Step &step);

	step::Time time;
	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	StructuralMechanics assembler;

	Matrix_Base<double> *K;
	Vector_Base<double> *U, *R, *f, *x, *dirichlet;

	SparseMatrixBuilder<double> *builder;
	LinearSystemSolver<double> *solver;

protected:
	bool checkDisplacement(step::Step &step);

	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_ */
