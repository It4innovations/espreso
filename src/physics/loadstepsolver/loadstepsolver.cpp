
#include "loadstepsolver.h"
#include "esinfo/stepinfo.h"
#include "esinfo/eslog.h"
#include "physics/substepsolver/substepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(LinearSystem *system, SubStepSolver *subStepSolver, double duration)
: system(system), subStepSolver(subStepSolver)
{
	step::time.start = step::time.current;
	step::time.final = step::time.current + duration;
	step::time.shift = duration;
}
