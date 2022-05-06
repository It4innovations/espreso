
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_

#include "harmonic.h"

namespace espreso {

struct HarmonicReal: public Harmonic {

	HarmonicReal(HarmonicSolverConfiguration &configuration, int dofs): Harmonic(configuration, dofs) {}

	void composeSystem(step::Frequency &frequency, LinearSystem<double> *system);
	void extractSolution(step::Frequency &frequency, LinearSystem<double> *system);
};

}



#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_ */
