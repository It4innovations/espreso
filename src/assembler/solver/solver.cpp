
#include "solver.h"

#include "../../basis/utilities/utils.h"
#include "../../configuration/environment.h"
#include "../../configuration/output.h"
#include "../../configuration/solver/espresooptions.h"
#include "../physics/physics.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../../mesh/structures/mesh.h"

#include "../../output/vtk/vtk.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

Solver::Solver(
		Mesh *mesh,
		std::vector<Physics*> &physics,
		std::vector<Instance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: physics(physics), instances(instances), linearSolvers(linearSolvers), _mesh(mesh), _store(store)
{
	_timeStatistics = new TimeEval("Solver Overall Timing");
	_timeStatistics->totalTime.startWithBarrier();
}

Solver::~Solver()
{
	delete _timeStatistics;
}


void Solver::storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description)
{
	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing " << description;
		for (size_t d = 0; d < matrices.size(); d++) {
			std::ofstream os(Logging::prepareFile(d, name));
			os << matrices[d];
			os.close();
		}
	}
}

void Solver::storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description)
{
	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing " << description;
		for (size_t d = 0; d < vectors.size(); d++) {
			std::ofstream os(Logging::prepareFile(d, name));
			os << vectors[d];
			os.close();
		}
	}
}

void Solver::storeSolution(const Step &step)
{
	for (size_t i = 0; i < instances.size(); i++) {
		for (size_t s = 0; s < instances[i]->solutions.size(); s++) {
			_store->storeValues(instances[i]->solutions[s]->name, instances[i]->solutions[i]->properties, instances[i]->solutions[s]->data, instances[i]->solutions[s]->eType);
		}
	}
}

void Solver::assembleStiffnessMatrices(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	std::for_each(physics.begin(), physics.end(), [&] (Physics *p) { p->assembleStiffnessMatrices(step); });
	timePhysics.endWithBarrier(); _timeStatistics->addEvent(timePhysics);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->K, instances.size() > 1 ? "K" + std::to_string(i) + "_" : "K", "stiffness matrices");
		storeData(step, instances[i]->f, instances.size() > 1 ? "f" + std::to_string(i) + "_" : "f", "RHS");
	}
}

void Solver::assembleResidualForces(const Step &step)
{
	ESINFO(PROGRESS2) << "Subtract residual forces.";
	TimeEvent timeSub("Subtract residual forces"); timeSub.start();
	for (size_t i = 0; i < instances.size(); i++) {
		physics[i]->assembleResidualForces(step);
	}
	timeSub.endWithBarrier(); _timeStatistics->addEvent(timeSub);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->f, instances.size() > 1 ? "f" + std::to_string(i) + "_" : "f", "RHS");
	}
}

void Solver::makeStiffnessMatricesRegular(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble matrices R1, R2 and RegMat.";
	TimeEvent timeReg("Make K regular"); timeReg.start();
	for (size_t i = 0; i < instances.size(); i++) {
		physics[i]->makeStiffnessMatricesRegular(linearSolvers[i]->configuration.regularization);
	}
	timeReg.endWithBarrier(); _timeStatistics->addEvent(timeReg);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->N1, instances.size() > 1 ? "R1" + std::to_string(i) + "_" : "R1", "R1");
		storeData(step, instances[i]->N2, instances.size() > 1 ? "R2" + std::to_string(i) + "_" : "R2", "R2");
		storeData(step, instances[i]->RegMat, instances.size() > 1 ? "RegMat" + std::to_string(i) + "_" : "RegMat", "regularization matrix");
	}
}

void Solver::assembleB1(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble B1.";
	TimeEvent timeConstrainsB1("Assemble B1"); timeConstrainsB1.startWithBarrier();
	for (size_t i = 0; i < instances.size(); i++) {
		for (size_t d = 0; d < instances[i]->domains; d++) {
			instances[i]->B1[d].cols = instances[i]->DOFs[d];
		}
		physics[i]->assembleB1(step, linearSolvers[i]->configuration.redundant_lagrange, linearSolvers[i]->configuration.scaling);
	}
	timeConstrainsB1.end(); _timeStatistics->addEvent(timeConstrainsB1);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->B1, instances.size() > 1 ? "B1" + std::to_string(i) + "_" : "B1", "Total FETI gluing matrices");
		storeData(step, instances[i]->B1duplicity, instances.size() > 1 ? "B1duplicity" + std::to_string(i) + "_" : "B1duplicity", "B1 duplicity");
		storeData(step, instances[i]->B1c, instances.size() > 1 ? "B1c" + std::to_string(i) + "_" : "B1c", "B1 c");
	}

	if (_store->configuration().gluing) {
		for (size_t i = 0; i < instances.size(); i++) {
			store::VTK::gluing(_store->configuration(), *_mesh, *instances[i], "B1", physics[i]->pointDOFs().size());
		}
	}
}

void Solver::subtractSolutionFromB1c(const Step &step)
{
	ESINFO(PROGRESS2) << "Subtract solution from B1c.";
	TimeEvent timesubtractB1c("Assemble B1"); timesubtractB1c.startWithBarrier();
	for (size_t i = 0; i < instances.size(); i++) {
		#pragma omp parallel for
		for (size_t d = 0; d < instances[i]->domains; d++) {
			for (size_t j = 0; j < instances[i]->B1[d].J_col_indices.size(); j++) {
				if (instances[i]->B1[d].I_row_indices[j] > instances[i]->block[Instance::CONSTRAINT::DIRICHLET]) {
					break;
				}
				instances[i]->B1c[d][j] -= instances[i]->primalSolution[d][instances[i]->B1[d].J_col_indices[j] - 1];
			}
		}
	}
	timesubtractB1c.end(); _timeStatistics->addEvent(timesubtractB1c);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->B1c, instances.size() > 1 ? "B1c" + std::to_string(i) + "_" : "B1c", "B1 c");
	}
}

void Solver::assembleB0(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble B0.";
	TimeEvent timeConstrainsB0("Assemble B0"); timeConstrainsB0.startWithBarrier();
	for (size_t i = 0; i < instances.size(); i++) {
		for (size_t d = 0; d < instances[i]->domains; d++) {
			instances[i]->B0[d].cols = instances[i]->DOFs[d];
		}

		if (linearSolvers[i]->configuration.method == ESPRESO_METHOD::TOTAL_FETI) {
			continue;
		}

		switch (linearSolvers[i]->configuration.B0_type) {
		case B0_TYPE::CORNERS:
			physics[i]->assembleB0FromCorners(step);
			break;
		case B0_TYPE::KERNELS:
			physics[i]->assembleB0FromKernels(step);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
	}
	timeConstrainsB0.end(); _timeStatistics->addEvent(timeConstrainsB0);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->B0, instances.size() > 1 ? "B0" + std::to_string(i) + "_" : "B0", "Hybrid Totoal FETI gluing matrices");
	}
}

void Solver::lineSearch(const std::vector<std::vector<double> > &previous, std::vector<std::vector<double> > &delta, Physics *physics, const Step &step)
{
	auto multiply = [] (const std::vector<std::vector<double> > &v1, const std::vector<std::vector<double> > &v2) {
		double cmul = 0, gmul;

		#pragma omp parallel for
		for (size_t d = 0; d < v1.size(); d++) {
			double dmul = 0;
			for (size_t i = 0; i < v1[d].size(); i++) {
				dmul += v1[d][i] * v2[d][i];
			}
			#pragma omp atomic
			cmul += dmul;
		}
		MPI_Allreduce(&cmul, &gmul, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return gmul;
	};

	double a = 0, b = 1, alpha = 1;
	double fa, fb, fx, faStart;

	std::vector<std::vector<double> > solution = delta;
	std::vector<std::vector<double> > R0 = physics->instance()->R;

	for (size_t i = 0; i < 6; i++) {
		sumVectors(solution, previous, delta, 1, alpha); // solution = U + alpha * deltaU

		// physics->instance()->R = K * solution
		solution.swap(physics->instance()->primalSolution);
		physics->assembleResidualForces(step);
		solution.swap(physics->instance()->primalSolution);

		if (i == 0) {
			faStart = multiply(delta, physics->instance()->f); // deltaU * RHS_i
			fb = -multiply(delta, physics->instance()->R); // deltaU * -(K * solution)
			if (faStart * fb > 0) {
				R0.swap(physics->instance()->R);
				return;
			}
			fa = faStart;
		} else {
			fx = -multiply(delta, physics->instance()->R); // deltaU * -(K * solution)
			if (fabs(fx) <= 0.5 * faStart) {
				alpha = a - fa * ((b - a ) / (fb - fa));
				break;
			}
			if (fa * fx < 0) {
				b = alpha;
				fb = fx;
			} else if (fb * fx < 0) {
				a = alpha;
				fa = fx;
			}
		}

		alpha = a - fa * ((b - a ) / (fb - fa));
	}

	if (alpha < 0.1) {
		alpha = 0.1;
	}
	if (alpha > 1) {
		alpha = 1;
	}

	sumVectors(solution, previous, delta, 0, alpha);
	solution.swap(delta);
}

void Solver::sumVectors(std::vector<std::vector<double> > &result, const std::vector<std::vector<double> > &a, const std::vector<std::vector<double> > &b, double alpha, double beta)
{
	#pragma omp parallel for
	for (size_t d = 0; d < a.size(); d++) {
		for (size_t i = 0; i < a[d].size(); i++) {
			result[d][i] = alpha * a[d][i] + beta * b[d][i];
		}
	}
}

void Solver::processSolution(const Step &step)
{
	TimeEvent store("Store solution"); store.startWithBarrier();
	for (size_t i = 0; i < instances.size(); i++) {
		physics[i]->processSolution(step);
	}
	store.end(); _timeStatistics->addEvent(store);

	for (size_t i = 0; i < instances.size(); i++) {
		storeData(step, instances[i]->primalSolution, instances.size() > 1 ? "solution" + std::to_string(i) + "_" : "solution", "solution");
	}
}

void Solver::initLinearSolver()
{
	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	for (size_t i = 0; i < linearSolvers.size(); i++) {
		linearSolvers[i]->steel(instances[i]);
		linearSolvers[i]->init(_mesh->neighbours());
	}
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::startLinearSolver()
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	for (size_t i = 0; i < instances.size(); i++) {
		linearSolvers[i]->Solve(instances[i]->f, instances[i]->primalSolution);
	}
	timeSolve.endWithBarrier(); _timeStatistics->addEvent(timeSolve);
}

void Solver::finalizeLinearSolver()
{
	std::for_each(linearSolvers.begin(), linearSolvers.end(), [&] (LinearSolver *s) { s->finilize(); });

	_store->finalize();

	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

