
#ifndef SOLVER_H_
#define SOLVER_H_

#include "../configuration.h"

template<class TAssembler>
class Solver {

public:
	Solver(const Instance &instance): _assembler(instance) { };

	void init();

	void solve( eslocal steps );

private:
	TAssembler _assembler;
};

#include "solver.hpp"


#endif /* SOLVER_H_ */
