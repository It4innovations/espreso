
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HarmonicSolverConfiguration;
namespace step { struct Frequency; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

struct Harmonic {

	Harmonic(HarmonicSolverConfiguration &configuration, int dofs);
	~Harmonic();

	void initFrequency(step::Frequency &frequency);
	void nextFrequency(step::Frequency &frequency);

	template <typename Assembler, typename Solver>
	void init(LinearSystem<Assembler, Solver> *system)
	{
		system->setMapping(K = system->assembler.A->copyPattern());
		system->setMapping(M = system->assembler.A->copyPattern());
		system->setMapping(C = system->assembler.A->copyPattern());
		system->setMapping(re.f = system->assembler.b->copyPattern());
		system->setMapping(im.f = system->assembler.b->copyPattern());
		system->setMapping(re.x = system->assembler.b->copyPattern());
		system->setMapping(im.x = system->assembler.b->copyPattern());
		system->setDirichletMapping(re.dirichlet = system->assembler.dirichlet->copyPattern());
		system->setDirichletMapping(im.dirichlet = system->assembler.dirichlet->copyPattern());
	}

	HarmonicSolverConfiguration &configuration;

	int dofs;
	Matrix_Base<double> *K, *M, *C;
	struct {
		Vector_Base<double> *f, *x, *dirichlet;
	} re, im;
};

}

#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_H_ */
