
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_

#include "harmonic.h"

namespace espreso {

struct AX_HarmonicReal: public AX_Harmonic {

	AX_HarmonicReal(HarmonicSolverConfiguration &configuration, int dofs): AX_Harmonic(configuration, dofs) {}

	void composeSystem(step::Frequency &frequency, AX_LinearSystem<double> *system);
	void extractSolution(step::Frequency &frequency, AX_LinearSystem<double> *system);
};

}



#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_ */
