
#include "solver.h"

template<class TAssembler>
void Solver<TAssembler>::solve()
{
	int time = 100;
	for (int i = 0; i < time; i++) {
		_assembler.update();
		_assembler.solve();
	}
}



