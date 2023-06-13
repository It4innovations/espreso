
#ifndef SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_
#define SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_

#include "analysis.h"
#include "analysis/linearsystem/linearsystem.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsSteadyStateNonLinear: public Analysis {

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

	LinearSystem<double> *system;

protected:
	bool checkDisplacement(step::Step &step);

	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_ */
