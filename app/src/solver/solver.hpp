
#include "solver.h"

template<class TAssembler>
void Solver<TAssembler>::solve(eslocal steps)
{
	_assembler.init();

	for (int i = 0; i < steps; i++) {
		_assembler.pre_solve_update();
		_assembler.solve();
		_assembler.post_solve_update();
	}

	_assembler.finalize();

}

template<class TAssembler>
void Solver<TAssembler>::init()
{

	//_assembler.init();

}

