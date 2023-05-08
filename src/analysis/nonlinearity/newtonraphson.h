
#ifndef SRC_ANALYSIS_NONLINEARITY_NEWTONRAPHSON_H_
#define SRC_ANALYSIS_NONLINEARITY_NEWTONRAPHSON_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; struct Time; }

template<typename T> struct Vector_Base;
template<typename T> struct Matrix_Base;

class NonLinearSolverConfiguration;
class HeatTransfer;
class StructuralMechanics;
class SteadyState;

class NewtonRaphson {

public:
	NewtonRaphson(NonLinearSolverConfiguration &configuration);
	~NewtonRaphson();

	void init(LinearSystem<double> *system);
	bool run(step::Step &step, step::Time &time, HeatTransfer &assembler, SteadyState &scheme, LinearSystem<double> *system);
	bool run(step::Step &step, step::Time &time, StructuralMechanics &assembler, SteadyState &scheme, LinearSystem<double> *system);

	NonLinearSolverConfiguration &configuration;
	Vector_Base<double> *U, *R;

protected:
	bool checkTemp(step::Step &step, HeatTransfer &assembler, SteadyState &scheme, LinearSystem<double> *system);
	bool checkDisplacement(step::Step &step, StructuralMechanics &assembler, SteadyState &scheme, LinearSystem<double> *system);
};

}


#endif /* SRC_ANALYSIS_NONLINEARITY_NEWTONRAPHSON_H_ */
