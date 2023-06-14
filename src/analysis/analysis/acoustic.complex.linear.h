
#ifndef SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_

#include "analysis.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/linearsystem/linearsystem.h"

#include <complex>

namespace espreso {

struct AcousticConfiguration;
struct AcousticLoadStepConfiguration;

class AcousticComplexLinear: public Analysis {

public:
	AcousticComplexLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration);

	void analyze();
	void run(step::Step &step);

	step::Frequency frequency;
	AcousticConfiguration &settings;
	AcousticLoadStepConfiguration &configuration;

	Acoustic assembler;

	Matrix_Base<double> *K, *M, *C;
	struct {
		Vector_Base<double> *f, *x, *dirichlet;
	} re, im;

	LinearSystem<double, std::complex<double> > *system;

protected:
	void storeSystem(step::Step &step);
	void storeSolution(step::Step &step);
};

}





#endif /* SRC_ANALYSIS_ANALYSIS_ACOUSTIC_COMPLEX_LINEAR_H_ */
