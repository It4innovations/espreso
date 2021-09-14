
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_

#include "analysis/linearsystem/linearsystem.h"

#include <complex>

namespace espreso {

struct HarmonicSolverConfiguration;
namespace step { struct Frequency; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

struct AX_HarmonicComplex {

	AX_HarmonicComplex(HarmonicSolverConfiguration &configuration, int dofs);
	~AX_HarmonicComplex();

	void initFrequency(step::Frequency &frequency);
	void nextFrequency(step::Frequency &frequency);

	void init(AX_LinearSystem<double, std::complex<double> > *system);

	void composeSystem(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system);
	void composeDirichlet(AX_LinearSystem<double, std::complex<double> > *system);

	void extractSolution(AX_LinearSystem<double, std::complex<double> > *system);

	void storeScheme(step::Frequency &frequency);
	void storeSolution(step::Frequency &frequency);

	HarmonicSolverConfiguration &configuration;

	int dofs;
	Matrix_Base<double> *K, *M, *C;
	struct {
		Vector_Base<double> *f, *x;
	} re, im;
	

	//Vector_Base<std::complex<double>> *f, *x;
};

}



#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_ */
