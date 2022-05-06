
#ifndef SRC_ANALYSIS_SCHEME_STEADYSTATE_H_
#define SRC_ANALYSIS_SCHEME_STEADYSTATE_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; struct Time; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

struct SteadyState {

	SteadyState();
	~SteadyState();

	void setTime(step::Time &time, double current);

	void init(LinearSystem<double> *system);

	void composeSystem(step::Step &step, LinearSystem<double> *system);
	void extractSolution(step::Step &step, LinearSystem<double> *system);

	Matrix_Base<double> *K;
	Vector_Base<double> *f, *x, *dirichlet;
};

}

#endif /* SRC_ANALYSIS_SCHEME_STEADYSTATE_H_ */
