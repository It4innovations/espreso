
#include "newtonraphson.h"

using namespace espreso;


void AX_NewtonRaphson::init(AX_LinearSystem<double> *system)
{
	if (this->system != nullptr) {
		// transform
	} else {
		K = system->solver.A->copyPattern();
		U = system->solver.b->copyPattern();
		R = system->solver.b->copyPattern();
		f = system->solver.b->copyPattern();
		BC = system->solver.b->copyPattern();
		f = system->solver.b->copyPattern();
		lsSolution = system->solver.b->copyPattern();
		lsRHS = system->solver.b->copyPattern();
		lsResidual = system->solver.b->copyPattern();
	}
	this->system = system;
}

