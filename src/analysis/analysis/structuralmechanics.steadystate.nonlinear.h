
#ifndef SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_
#define SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_

#include "analysis.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"
#include "analysis/assembler/module/structuralmechanics.h"
#include "analysis/nonlinearity/newtonraphson.h"

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
	NewtonRaphson solver;
	SteadyState scheme;

	LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_ */
