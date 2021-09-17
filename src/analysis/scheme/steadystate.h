
#ifndef SRC_ANALYSIS_SCHEME_STEADYSTATE_H_
#define SRC_ANALYSIS_SCHEME_STEADYSTATE_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; struct Time; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

struct AX_SteadyState {

	AX_SteadyState();
	~AX_SteadyState();

	void setTime(step::Time &time, double current);

	void init(AX_LinearSystem<double> *system);

	void composeSystem(step::Step &step, AX_LinearSystem<double> *system);
	void extractSolution(step::Step &step, AX_LinearSystem<double> *system);

	Matrix_Base<double> *K;
	Vector_Base<double> *f, *x, *dirichlet;
};

}

#endif /* SRC_ANALYSIS_SCHEME_STEADYSTATE_H_ */
