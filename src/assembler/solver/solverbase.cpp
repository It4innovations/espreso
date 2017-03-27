
#include "solverbase.h"

#include "../../basis/logging/timeeval.h"

using namespace espreso;

SolverBase::SolverBase(const std::string &name, const std::string &physicsName, Mesh *mesh)
: _name(name), _mesh(mesh)
{
	_timeStatistics = new TimeEval(physicsName + " solved by " + _name + " solver overall timing");
	_timeStatistics->totalTime.startWithBarrier();
}

SolverBase::~SolverBase()
{
	delete _timeStatistics;
}


