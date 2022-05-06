
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_

#include "harmonic.h"

namespace espreso {

struct HarmonicComplex: public Harmonic {

	HarmonicComplex(HarmonicSolverConfiguration &configuration, int dofs): Harmonic(configuration, dofs) {}

	void composeSystem(step::Frequency &frequency, LinearSystem<double, std::complex<double> > *system);
	void extractSolution(step::Frequency &frequency, LinearSystem<double, std::complex<double> > *system);
};

}

#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_ */
