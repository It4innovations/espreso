
#include "linear.h"

#include "analysis/linearsolver/directsolver.h"
#include "analysis/linearsolver/fetisolver.h"
#include "analysis/linearsolver/multigridsolver.h"

using namespace espreso;

void AX_Linear::init(DirectSolver<double> *solver)
{
	printf("init Linear solver for Direct<double>\n");
	this->solver = solver;
}

void AX_Linear::init(FETISolver<double> *solver)
{
	printf("init Linear solver for FETI<double>\n");
	this->solver = solver;
}

void AX_Linear::init(MultigridSolver<double> *solver)
{
	printf("init Linear solver for Multigrid<double>\n");
	this->solver = solver;
}
