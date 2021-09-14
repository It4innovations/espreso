
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HarmonicSolverConfiguration;
namespace step { struct Frequency; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

struct AX_HarmonicReal {

	AX_HarmonicReal(HarmonicSolverConfiguration &configuration, int dofs);
	~AX_HarmonicReal();

	void initFrequency(step::Frequency &frequency);
	void nextFrequency(step::Frequency &frequency);

	void init(AX_LinearSystem<double> *system);

	void composeSystem(step::Frequency &frequency, AX_LinearSystem<double> *system);
	void composeDirichlet(AX_LinearSystem<double> *system);

	void extractSolution(AX_LinearSystem<double> *system);

	void storeScheme(step::Frequency &frequency);
	void storeSolution(step::Frequency &frequency);

	HarmonicSolverConfiguration &configuration;

	int dofs;
	Matrix_Base<double> *K, *M, *C;
	struct {
		Vector_Base<double> *f, *x;
	} re, im;
};

}



#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_ */
