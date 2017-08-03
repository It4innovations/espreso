
#include "linear.h"

#include "../assembler.h"
#include "../loadstep/loadstepsolver.h"
#include "../../instance.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Assembler &assembler)
: TimeStepSolver("LINEAR", assembler)
{

}

void LinearTimeStep::solve(Step &step, LoadStepSolver &loadStepSolver)
{
	_assembler.solve(step, loadStepSolver.updateStructuralMatrices(step, Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
}




