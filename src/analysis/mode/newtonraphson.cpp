
#include "newtonraphson.h"

#include "analysis/linearsolver/directsolver.h"
#include "analysis/linearsolver/fetisolver.h"
#include "analysis/linearsolver/multigridsolver.h"

using namespace espreso;

void AX_NewtonRaphson::init(DirectSolver<double> *solver)
{
	if (this->solver != nullptr) {

	} else {

	}
	this->solver = solver;
}

void AX_NewtonRaphson::init(FETISolver<double> *solver)
{
	if (this->solver != nullptr) {
		// transform
	} else {
		K = solver->getMatrixPattern();
		U = solver->getVectorPattern();
		R = solver->getVectorPattern();
		f = solver->getVectorPattern();
		BC = solver->getVectorPattern();
		f = solver->getVectorPattern();
		lsSolution = solver->getVectorPattern();
		lsRHS = solver->getVectorPattern();
		lsResidual = solver->getVectorPattern();
	}
	this->solver = solver;
}

void AX_NewtonRaphson::init(MultigridSolver<double> *solver)
{
	if (this->solver != nullptr) {

	} else {

	}
	this->solver = solver;
}
