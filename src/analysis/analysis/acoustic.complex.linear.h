
#ifndef SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/scheme/harmonic.complex.h"
#include "analysis/linearsystem/linearsystem.h"

#include <complex>

namespace espreso {

struct AcousticGlobalSettings;
struct AcousticLoadStepConfiguration;

class AX_AcousticComplexLinear: public Analysis {

public:
	AX_AcousticComplexLinear(AcousticGlobalSettings &gsettings, AcousticLoadStepConfiguration &configuration);

	void init();
	void run(step::Step &step);

	AcousticGlobalSettings &gsettings;
	AcousticLoadStepConfiguration &configuration;

	AX_Acoustic assembler;
	AX_HarmonicComplex scheme;

	AX_LinearSystem<std::complex<double> > *system;
};

}





#endif /* SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_ */
