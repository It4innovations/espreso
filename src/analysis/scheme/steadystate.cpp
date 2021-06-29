
#include "steadystate.h"

#include "analysis/linearsolver/directsolver.h"
#include "analysis/linearsolver/fetisolver.h"
#include "analysis/linearsolver/multigridsolver.h"

using namespace espreso;

void AX_SteadyState::init(DirectSolver<double> *solver)
{

}

void AX_SteadyState::init(FETISolver<double> *solver)
{
	printf("init STEADY STATE scheme for FETI<double>\n");
	K = solver->getMatrixPattern();
	f = solver->getVectorPattern();
}

void AX_SteadyState::init(MultigridSolver<double> *solver)
{

}
