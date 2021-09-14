
#ifndef SRC_ANALYSIS_SOLVER_NEWTONRAPHSON_H_
#define SRC_ANALYSIS_SOLVER_NEWTONRAPHSON_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; struct Time; }

template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

class NonLinearSolverConfiguration;
class AX_HeatTransfer;
class AX_SteadyState;

class AX_NewtonRaphson {

public:
	AX_NewtonRaphson(NonLinearSolverConfiguration &configuration);
	~AX_NewtonRaphson();

	void init(AX_LinearSystem<double> *system);
	bool run(step::Step &step, step::Time &time, AX_HeatTransfer &assembler, AX_SteadyState &scheme, AX_LinearSystem<double> *system);

protected:
	bool checkTemp(step::Step &step, AX_HeatTransfer &assembler, AX_SteadyState &scheme, AX_LinearSystem<double> *system);

	NonLinearSolverConfiguration &configuration;

	Vector_Base<double> *U, *R;
};

}


#endif /* SRC_ANALYSIS_SOLVER_NEWTONRAPHSON_H_ */
