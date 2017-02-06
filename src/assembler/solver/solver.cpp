
#include "solver.h"

#include "../../basis/utilities/utils.h"
#include "../../config/environment.h"
#include "../../config/solverespresooptions.h"
#include "../../config/output.h"

#include "../physics/physics.h"
#include "../step.h"
#include "../instance.h"

#include "../../mesh/structures/mesh.h"

#include "../../output/vtk/vtk.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

Solver::Solver(
		Mesh *mesh,
		std::vector<NewPhysics*> &physics,
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

void Solver::assembleStiffnessMatrices(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	std::for_each(physics.begin(), physics.end(), [&] (NewPhysics *p) { p->assembleStiffnessMatrices(step); });
	timePhysics.endWithBarrier(); _timeStatistics->addEvent(timePhysics);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing stiffness matrices and RHS";
		for (size_t i = 0; i < instances.size(); i++) {
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "K" + std::to_string(i) + "_" : "K").c_str());
				os << instances[i]->K[d];
				os.close();
			}
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "f" + std::to_string(i) + "_" : "f").c_str());
				os << instances[i]->f[d];
				os.close();
			}
		}
	}
}

void Solver::makeStiffnessMatricesRegular()
{
	ESINFO(PROGRESS2) << "Assemble matrices R1, R2 and RegMat.";
	TimeEvent timeReg("Make K regular"); timeReg.start();
	for (size_t i = 0; i < instances.size(); i++) {
		physics[i]->makeStiffnessMatricesRegular(linearSolvers[i]->configuration.regularization);
	}
	timeReg.endWithBarrier(); _timeStatistics->addEvent(timeReg);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing R1, R2, RegMat";
		for (size_t i = 0; i < instances.size(); i++) {
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "R1" + std::to_string(i) + "_" : "R1").c_str());
				os << instances[i]->R1[d];
				os.close();
			}
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "R2" + std::to_string(i) + "_" : "R2").c_str());
				os << instances[i]->R2[d];
				os.close();
			}
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "RegMat" + std::to_string(i) + "_" : "RegMat").c_str());
				os << instances[i]->RegMat[d];
				os.close();
			}
		}
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

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing B1";
		for (size_t i = 0; i < instances.size(); i++) {
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "B1" + std::to_string(i) + "_" : "B1").c_str());
				os << instances[i]->B1[d];
				os.close();
			}
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "B1duplicity" + std::to_string(i) + "_" : "B1duplicity").c_str());
				os << instances[i]->B1duplicity[d];
				os.close();
			}
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "B1c" + std::to_string(i) + "_" : "B1c").c_str());
				os << instances[i]->B1c[d];
				os.close();
			}
		}
	}
	if (_store->configuration().gluing) {
		for (size_t i = 0; i < instances.size(); i++) {
			store::VTK::gluing(_store->configuration(), *_mesh, *instances[i], "B1", physics[i]->pointDOFs().size());
		}
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

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing B0";
		for (size_t i = 0; i < instances.size(); i++) {
			for (size_t d = 0; d < instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, instances.size() > 1 ? "B0" + std::to_string(i) + "_" : "B0").c_str());
				os << instances[i]->B0[d];
				os.close();
			}
		}
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

void Solver::startLinearSolver(const Step &step)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	for (size_t i = 0; i < instances.size(); i++) {
		linearSolvers[i]->Solve(instances[i]->f, instances[i]->primalSolution);
		physics[i]->storeSolution(step, instances[i]->primalSolution, _store);
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

