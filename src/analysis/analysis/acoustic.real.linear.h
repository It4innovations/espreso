
#ifndef SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/scheme/harmonic.real.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct AcousticGlobalSettings;
struct AcousticLoadStepConfiguration;

class AX_AcousticRealLinear: public Analysis {

public:
	AX_AcousticRealLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration);

	void analyze();
	void run(step::Step &step);

	AcousticConfiguration &settings;
	AcousticLoadStepConfiguration &configuration;

	AX_Acoustic assembler;
	AX_HarmonicReal scheme;

	AX_LinearSystem<double> *system;
};

}




#endif /* SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_ */
