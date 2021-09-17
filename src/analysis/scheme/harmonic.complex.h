
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_

#include "harmonic.h"

namespace espreso {

struct AX_HarmonicComplex: public AX_Harmonic {

	AX_HarmonicComplex(HarmonicSolverConfiguration &configuration, int dofs): AX_Harmonic(configuration, dofs) {}

	void composeSystem(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system);
	void extractSolution(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system);
};

}

#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_ */
