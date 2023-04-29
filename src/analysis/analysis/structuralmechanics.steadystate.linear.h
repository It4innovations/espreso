
#ifndef SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_

#include "analysis.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsSteadyStateLinear: public Analysis {

public:
	StructuralMechanicsSteadyStateLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanicsSteadyStateLinear();

	void analyze();
	void run(step::Step &step);

	step::Time time;
	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	StructuralMechanics assembler;
	SteadyState scheme;

	LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_ELASTICITY_STEADYSTATE_LINEAR_H_ */
