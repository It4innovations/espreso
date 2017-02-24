
#include "transientfirstorderimplicit.h"

#include "../step.h"
#include "../instance.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		store::ResultStore* store)
: Solver("TRANSIENT FIRST ORDER IMPLICIT", mesh, physics, linearSolver, store, Matrices::NONE)
{

}

void TransientFirstOrderImplicit::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	assembleMatrices(step, Matrices::M | Matrices::K | Matrices::f);
	composeGluing(step, Matrices::B1);

	// TODO: B0 without kernels
	composeGluing(step, Matrices::B0);

	// case
	// Crank Nicolson alpha = 0.5
	// Galerkin  alpha = 2 / 3
	// Backward difference alpha = 1

	// K += 1 / alpha deltaT * M

	// f += M * ( 1 / alpha deltaT * prevSolution + (1 - alpha) / alpha * prevRychlost)

	// solve()

	// rychlost = (1 / alpha deltaT) * (solution - prevSolution) - (1 - alpha) / alpha * prevRychlost
}

void TransientFirstOrderImplicit::init(Step &step)
{

}

void TransientFirstOrderImplicit::preprocess(Step &step)
{
}

void TransientFirstOrderImplicit::solve(Step &step)
{

}

void TransientFirstOrderImplicit::postprocess(Step &step)
{
}

void TransientFirstOrderImplicit::finalize(Step &step)
{
}




