
#include "transientfirstorderimplicit.h"


#include "../../configuration/physics/transientsolver.h"

#include "../step.h"
#include "../instance.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		store::ResultStore* store,
		const TransientSolver &configuration)
: Solver("TRANSIENT FIRST ORDER IMPLICIT", mesh, physics, linearSolver, store, Matrices::NONE),
  _configuration(configuration)
{

}

void TransientFirstOrderImplicit::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	double alpha = 0;
	switch (_configuration.method) {
	case TransientSolver::METHOD::CRANK_NICOLSON:
		alpha = 0.5;
		break;
	case TransientSolver::METHOD::GALERKIN:
		alpha = 2 / 3;
		break;
	case TransientSolver::METHOD::BACKWARD_DIFF:
		alpha = 1;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not supported first order implicit solver method.";
	}

	std::vector<std::vector<double> > u(instance->domains), deltaU(instance->domains);
	std::vector<std::vector<double> > v(instance->domains);
	std::vector<std::vector<double> > x(instance->domains), y(instance->domains);

	for (size_t d = 0; d < instance->domains; d++) {
		u[d].resize(instance->DOFs[d]);
		deltaU[d].resize(instance->DOFs[d]);
		v[d].resize(instance->DOFs[d]);
		x[d].resize(instance->DOFs[d]);
		y[d].resize(instance->DOFs[d]);
	}

	for (step.iteration = 0; step.iteration < 1; step.iteration++) {
		assembleMatrices(step, Matrices::K | Matrices::M | Matrices::f);
		composeGluing(step, Matrices::B1 | Matrices::B0); // TODO: B0 without kernels

		sum(step, Matrices::K, Matrices::M, 1, 1 / (alpha * _configuration.time_step));
		sum(x, u, v, 1 / (alpha * _configuration.time_step), (1 - alpha) / alpha, "x", "u", "v");
		multiply(step, Matrices::M, x, y, 1, "x", "y");
		sum(step, Matrices::f, y, 1, 1, "y");

		if (step.iteration) {
			updateLinearSolver(Matrices::K | Matrices::M | Matrices::f | Matrices::B1c);
		} else {
			initLinearSolver();
		}
		runLinearSolver();

		sum(deltaU, physics->instance()->primalSolution, u, 1, -1);
		sum(v, v, deltaU, 1 / (alpha * _configuration.time_step), (1 - alpha) / alpha, "v", "v", "deltaU");
		sum(step, Matrices::primar, deltaU, 1, 1, "deltaU");

		processSolution(step);
		storeSolution(step);
	}
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




