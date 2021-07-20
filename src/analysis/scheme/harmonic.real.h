
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_

namespace espreso {

struct HarmonicSolverConfiguration;
namespace step { struct Frequency; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;
template<typename T> struct AX_LinearSystem;

struct AX_HarmonicReal {

	AX_HarmonicReal(HarmonicSolverConfiguration &configuration, int dofs)
	: configuration(configuration), dofs(dofs), K{}, M{}, C{}, re{}, im{} {}

	void initFrequency(step::Frequency &frequency);
	void nextFrequency(step::Frequency &frequency);

	void init(AX_LinearSystem<double> *system);

	void composeSystem(step::Frequency &frequency, AX_LinearSystem<double> *system);
	void composeDirichlet(AX_LinearSystem<double> *system);

	void extractSolution(AX_LinearSystem<double> *system);

	void storeScheme(step::Frequency &frequency);
	void storeSolution(step::Frequency &frequency);

	HarmonicSolverConfiguration &configuration;

	struct Fragment {
		Vector_Base<double> *f, *x;
	};

	int dofs;
	Matrix_Base<double> *K, *M, *C;
	Fragment re, im;
};

}



#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_REAL_H_ */
