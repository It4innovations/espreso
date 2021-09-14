
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_

#include <complex>

namespace espreso {

struct HarmonicSolverConfiguration;
namespace step { struct Frequency; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;
template<typename T> struct AX_LinearSystem;

using AX_ComplexLinearSystem = AX_LinearSystem<std::complex<double> >;

struct AX_HarmonicComplex {

	AX_HarmonicComplex(HarmonicSolverConfiguration &configuration, int dofs);
	~AX_HarmonicComplex();

	void initFrequency(step::Frequency &frequency);
	void nextFrequency(step::Frequency &frequency);

	void init(AX_ComplexLinearSystem *system);

	void composeSystem(step::Frequency &frequency, AX_ComplexLinearSystem *system);
	void composeDirichlet(AX_ComplexLinearSystem *system);

	void extractSolution(AX_ComplexLinearSystem *system);

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
