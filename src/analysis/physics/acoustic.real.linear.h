
#ifndef SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/acoustic.h"
//#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct AcousticGlobalSettings;
struct AcousticLoadStepConfiguration;

class AcousticRealLinear: public Physics {

public:
	AcousticRealLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration);

	void analyze(step::Step &step);
	void run(step::Step &step);

	step::Frequency frequency;
	AcousticConfiguration &settings;
	AcousticLoadStepConfiguration &configuration;

	Acoustic assembler;

	Matrix_Base<double> *K, *M, *C;
	struct {
		Vector_Base<double> *f, *x, *dirichlet;
	} re, im;

//	LinearSystem<double> *system;

protected:
	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}




#endif /* SRC_ANALYSIS_ANALYSIS_ACOUSTIC_REAL_LINEAR_H_ */
