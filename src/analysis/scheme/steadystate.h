
#ifndef SRC_ANALYSIS_SCHEME_STEADYSTATE_H_
#define SRC_ANALYSIS_SCHEME_STEADYSTATE_H_

namespace espreso {

namespace step { struct Time; }
template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;
template<typename T> struct AX_LinearSystem;

struct AX_SteadyState {

	AX_SteadyState();
	~AX_SteadyState();

	void setTime(step::Time &time, double current);

	void init(AX_LinearSystem<double> *system);

	void composeSystem(AX_LinearSystem<double> *system);
	void extractSolution(AX_LinearSystem<double> *system);

	void storeScheme(step::Time &time);
	void storeSolution(step::Time &time);

	Matrix_Base<double> *K;
	Vector_Base<double> *f, *x;
};

}

#endif /* SRC_ANALYSIS_SCHEME_STEADYSTATE_H_ */
