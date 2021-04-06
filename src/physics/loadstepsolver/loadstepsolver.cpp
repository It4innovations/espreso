
#include "loadstepsolver.h"
#include "esinfo/stepinfo.h"
#include "esinfo/eslog.h"
#include "physics/substepsolver/substepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(LinearSystem *system, SubStepSolver *subStepSolver, double duration)
: _system(system), _subStepSolver(subStepSolver)
{
	step::time::start = step::time::current;
	step::time::final = step::time::current + duration;
	step::time::shift = duration;
}

void LoadStepSolver::run()
{
	step::substep = step::duplicate::offset;
	step::iteration = 0;

	while (!step::isLast()) {
		runNextSubstep();
		step::substep++;
		step::iteration = 0;
	}
}
