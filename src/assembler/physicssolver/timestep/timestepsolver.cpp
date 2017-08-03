
#include "timestepsolver.h"

using namespace espreso;

TimeStepSolver::TimeStepSolver(const std::string &description, Assembler &assembler)
: _description(description), _assembler(assembler)
{

}

std::string TimeStepSolver::description() const
{
	return _description;
}
