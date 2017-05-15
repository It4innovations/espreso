
#include "transientfirstorderimplicit.h"


#include "../../configuration/physics/transientsolver.h"

#include "../step.h"
#include "../instance.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/constants.h"

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		output::Store* store,
		const TransientSolver &configuration,
		double duration)
: Solver("TRANSIENT FIRST ORDER IMPLICIT", mesh, physics, linearSolver, store, duration, Matrices::NONE),
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
	case TransientSolver::METHOD::USER:
		alpha = _configuration.alpha;
		if (alpha <= 0 || alpha > 1) {
			ESINFO(GLOBAL_ERROR) << "Alpha has to be from interval (0, 1>.";
		}
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not supported first order implicit solver method.";
	}

	std::vector<std::vector<double> > u(instance->domains), deltaU(instance->domains);
	std::vector<std::vector<double> > v(instance->domains);
	std::vector<std::vector<double> > x(instance->domains), y(instance->domains);

	for (size_t d = 0; d < instance->domains; d++) {
		u[d].resize(instance->domainDOFCount[d]);
		deltaU[d].resize(instance->domainDOFCount[d]);
		v[d].resize(instance->domainDOFCount[d]);
		x[d].resize(instance->domainDOFCount[d]);
		y[d].resize(instance->domainDOFCount[d]);
	}

	double startTime = step.currentTime;
	step.timeStep = _configuration.time_step;
	step.currentTime += step.timeStep;

	for (step.substep = 0; step.currentTime < startTime + _duration; step.substep++) {
		ESINFO(PROGRESS2) << _name << " iteration " << step.substep + 1;
		if (step.substep) {
			updateMatrices(step, Matrices::K | Matrices::M | Matrices::f, instance->solutions);
		} else {
			assembleMatrices(step, Matrices::K | Matrices::M | Matrices::f);
		}
		composeGluing(step, Matrices::B1 | Matrices::B0); // TODO: B0 without kernels

		sum(instance->K, 1 / (alpha * _configuration.time_step), instance->M, "K += (1 / " + ASCII::alpha + ASCII::DELTA + "t) * M");
		sum(x, 1 / (alpha * _configuration.time_step), u, (1 - alpha) / alpha, v, "x = (1 / " + ASCII::alpha + ASCII::DELTA + ") * u + (1 - " + ASCII::alpha + ") / " + ASCII::alpha + " * v");

		multiply(y, instance->M, x, "y = M * x");
		sum(instance->f, 1, instance->f, 1, y, "f += y");

		if (step.substep) {
			updateLinearSolver(step, Matrices::K | Matrices::M | Matrices::f | Matrices::B1c);
		} else {
			initLinearSolver(step);
		}
		runLinearSolver(step);
		sum(deltaU, 1, instance->primalSolution, -1, u, ASCII::DELTA + "u = u_i - u_i_1");
		sum(v, 1 / (alpha * _configuration.time_step), deltaU, - (1 - alpha) / alpha, v, "v = (1 / " + ASCII::alpha + ASCII::DELTA + "t) * " + ASCII::DELTA + "u - (1 - " + ASCII::alpha + ") / " + ASCII::alpha + " * v");
		u = instance->primalSolution;

		processSolution(step);
		storeSolution(step);
		step.currentTime += step.timeStep;
	}
	finalizeLinearSolver(step);
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




