
#ifndef SOLVER_H_
#define SOLVER_H_

#include "../configuration.h"

template<class TAssembler>
class Solver {

public:
	Solver(const Instance &instance)
		: _configuration(instance.configuration()), _assembler(instance) { };

	void init();

	void solve( eslocal steps );

private:
	Configuration _configuration;
	TAssembler _assembler;
};

#include "solver.hpp"


#endif /* SOLVER_H_ */
