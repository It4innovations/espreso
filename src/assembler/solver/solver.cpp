
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
		std::vector<NewInstance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: _mesh(mesh), _physics(physics), _instances(instances), _linearSolvers(linearSolvers), _store(store)
{
	_timeStatistics = new TimeEval("Solver Overall Timing");
	_timeStatistics->totalTime.startWithBarrier();
}

Solver::~Solver()
{

}


void Solver::meshPreprocessing()
{
	ESINFO(PROGRESS2) << "Prepare mesh structures.";
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	for (size_t i = 0; i < _instances.size(); i++) {

		switch (_linearSolvers[i]->configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			_physics[i]->prepareTotalFETI();
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			switch (_linearSolvers[i]->configuration.B0_type) {
			case B0_TYPE::CORNERS:
				_physics[i]->prepareHybridTotalFETIWithCorners();
				break;
			case B0_TYPE::KERNELS:
				_physics[i]->prepareHybridTotalFETIWithKernels();
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown FETI method";
		}

	}
	timePreparation.endWithBarrier(); _timeStatistics->addEvent(timePreparation);

	TimeEvent timeStoring("Store mesh data"); timeStoring.start();
	_store->storeGeometry();
	timeStoring.endWithBarrier(); _timeStatistics->addEvent(timeStoring);
}

void Solver::assembleStiffnessMatrices(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble matrices K and RHS.";
	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	std::for_each(_physics.begin(), _physics.end(), [&] (NewPhysics *p) { p->assembleStiffnessMatrices(step); });
	timePhysics.endWithBarrier(); _timeStatistics->addEvent(timePhysics);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing stiffness matrices and RHS";
		for (size_t i = 0; i < _instances.size(); i++) {
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "K" + std::to_string(i) + "_" : "K").c_str());
				os << _instances[i]->K[d];
				os.close();
			}
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "f" + std::to_string(i) + "_" : "f").c_str());
				os << _instances[i]->f[d];
				os.close();
			}
		}
	}
}

void Solver::makeStiffnessMatricesRegular()
{
	ESINFO(PROGRESS2) << "Assemble matrices R1, R2 and RegMat.";
	TimeEvent timeReg("Make K regular"); timeReg.start();
	for (size_t i = 0; i < _instances.size(); i++) {
		_physics[i]->makeStiffnessMatricesRegular(_linearSolvers[i]->configuration.regularization);
	}
	timeReg.endWithBarrier(); _timeStatistics->addEvent(timeReg);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing R1, R2, RegMat";
		for (size_t i = 0; i < _instances.size(); i++) {
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "R1" + std::to_string(i) + "_" : "R1").c_str());
				os << _instances[i]->R1[d];
				os.close();
			}
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "R2" + std::to_string(i) + "_" : "R2").c_str());
				os << _instances[i]->R2[d];
				os.close();
			}
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "RegMat" + std::to_string(i) + "_" : "RegMat").c_str());
				os << _instances[i]->RegMat[d];
				os.close();
			}
		}
	}
}

void Solver::assembleB1(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble B1.";
	TimeEvent timeConstrainsB1("Assemble B1"); timeConstrainsB1.startWithBarrier();
	for (size_t i = 0; i < _instances.size(); i++) {
		for (size_t d = 0; d < _instances[i]->domains; d++) {
			_instances[i]->B1[d].cols = _instances[i]->DOFs[d];
		}
		_physics[i]->assembleB1(step, _linearSolvers[i]->configuration.redundant_lagrange, _linearSolvers[i]->configuration.scaling);
	}
	timeConstrainsB1.end(); _timeStatistics->addEvent(timeConstrainsB1);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing B1";
		for (size_t i = 0; i < _instances.size(); i++) {
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "B1" + std::to_string(i) + "_" : "B1").c_str());
				os << _instances[i]->B1[d];
				os.close();
			}
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "B1duplicity" + std::to_string(i) + "_" : "B1duplicity").c_str());
				os << _instances[i]->B1duplicity[d];
				os.close();
			}
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "B1c" + std::to_string(i) + "_" : "B1c").c_str());
				os << _instances[i]->B1c[d];
				os.close();
			}
		}
	}
	if (_store->configuration().gluing) {
		for (size_t i = 0; i < _instances.size(); i++) {
			store::VTK::gluing(_store->configuration(), *_mesh, *_instances[i], "B1", _physics[i]->pointDOFs().size());
		}
	}
}

void Solver::assembleB0(const Step &step)
{
	ESINFO(PROGRESS2) << "Assemble B0.";
	TimeEvent timeConstrainsB0("Assemble B0"); timeConstrainsB0.startWithBarrier();
	for (size_t i = 0; i < _instances.size(); i++) {
		for (size_t d = 0; d < _instances[i]->domains; d++) {
			_instances[i]->B0[d].cols = _instances[i]->DOFs[d];
		}

		if (_linearSolvers[i]->configuration.method == ESPRESO_METHOD::TOTAL_FETI) {
			continue;
		}

		switch (_linearSolvers[i]->configuration.B0_type) {
		case B0_TYPE::CORNERS:
			_physics[i]->assembleB0FromCorners(step);
			break;
		case B0_TYPE::KERNELS:
			_physics[i]->assembleB0FromKernels(step);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
	}
	timeConstrainsB0.end(); _timeStatistics->addEvent(timeConstrainsB0);

	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing B0";
		for (size_t i = 0; i < _instances.size(); i++) {
			for (size_t d = 0; d < _instances[i]->domains; d++) {
				std::ofstream os(Logging::prepareFile(d, _instances.size() > 1 ? "B0" + std::to_string(i) + "_" : "B0").c_str());
				os << _instances[i]->B0[d];
				os.close();
			}
		}
	}
}

void Solver::initLinearSolver()
{
	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	for (size_t i = 0; i < _linearSolvers.size(); i++) {
		_linearSolvers[i]->steel(_instances[i]);
		_linearSolvers[i]->init(_mesh->neighbours());
	}
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::startLinearSolver(const Step &step)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	for (size_t i = 0; i < _instances.size(); i++) {
		_linearSolvers[i]->Solve(_instances[i]->f, _instances[i]->primalSolution);
		_physics[i]->storeSolution(step, _instances[i]->primalSolution, _store);
	}
	timeSolve.endWithBarrier(); _timeStatistics->addEvent(timeSolve);
}

void Solver::finalizeLinearSolver()
{
	std::for_each(_linearSolvers.begin(), _linearSolvers.end(), [&] (LinearSolver *s) { s->finilize(); });

	_store->finalize();

	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

