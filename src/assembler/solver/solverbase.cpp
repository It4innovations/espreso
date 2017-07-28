
#include "solverbase.h"

#include "../../basis/logging/timeeval.h"
#include <iostream>

using namespace espreso;

SolverBase::SolverBase(const std::string &name, const std::string &physicsName, Mesh *mesh, double duration)
: _name(name), _mesh(mesh), _duration(duration)
{
	_timeStatistics = new TimeEval(physicsName + " solved by " + _name + " solver overall timing");
	_timeStatistics->totalTime.startWithBarrier();
}

SolverBase::~SolverBase()
{
	delete _timeStatistics;
}


