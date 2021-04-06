
#include "physicalsolver.h"
#include "esinfo/stepinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"

#include "loadstepsolver/pseudotimesteppingsolver.h"
#include "loadstepsolver/steadystatesolver.h"
#include "loadstepsolver/transientfirstorderimplicitsolver.h"
#include "loadstepsolver/transientsecondorderimplicitsolver.h"
#include "loadstepsolver/harmonicsolver.h"
#include "loadstepsolver/topologyoptimization.h"

#include "substepsolver/linearsubstepsolver.h"
#include "substepsolver/newtonraphsonsolver.h"

#include "system/fetisystem.h"
#include "system/hypresystem.h"
#include "system/mklpdsssystem.h"
#include "system/pardisosystem.h"
#include "system/superlusystem.h"
#include "system/wsmpsystem.h"
#include "system/builder/timebuilder.h"
#include "system/builder/harmonicbuilder.h"
#include "composer/distributed/nodes.uniform.distributed.composer.h"
#include "composer/distributed/faces.edges.uniform.distributed.composer.h"
#include "composer/feti/nodes.uniform.feti.composer.h"
#include "kernels/heattransfer2d.kernel.h"
#include "kernels/heattransfer3d.kernel.h"
#include "kernels/structuralmechanics2d.kernel.h"
#include "kernels/structuralmechanics3d.elasticity.kernel.h"
#include "kernels/structuralmechanics3d.harmonic.kernel.h"
#include "kernels/structuralmechanics3d.tdnns.kernel.h"
#include "kernels/utils/morphing.h"
#include "output/resultstore.h"

using namespace espreso;

static LinearSystem* getSystem(LinearSystem *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	LinearSystem *current = NULL;
	Kernel *kernel = NULL, *pKernel = previous ? previous->assembler()->composer->kernel : NULL;
	switch (dimension) {
	case DIMENSION::D2: kernel = new HeatTransfer2DKernel(dynamic_cast<HeatTransfer2DKernel*>(pKernel), physics, gsettings, loadStep); break;
	case DIMENSION::D3: kernel = new HeatTransfer3DKernel(dynamic_cast<HeatTransfer3DKernel*>(pKernel), physics, gsettings, loadStep); break;
	default: break;
	}

	switch (loadStep.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI: {
		FETISystem *system = new FETISystem(1, 1, loadStep.feti);
		system->assembler()->composer = new NodesUniformFETIComposer(loadStep.feti, kernel, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE: {
		HYPRESystem *system = new HYPRESystem(1, 1, loadStep.hypre);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: {
		MKLPDSSSystem *system = new MKLPDSSSystem(1, 1, loadStep.mklpdss);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: {
		PARDISOSystem *system = new PARDISOSystem(1, 1, loadStep.pardiso);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: {
		SuperLUSystem *system = new SuperLUSystem(1, 1, loadStep.superlu);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::WSMP: {
		WSMPSystem *system = new WSMPSystem(1, 1, loadStep.wsmp);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), 1);
		current = system;
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unknown linear solver.\n");
	}

	current->builder = new TimeBuilder();

	current->init();
	return current;
}

static LinearSystem* getSystem(LinearSystem *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	auto tdnns = [&] () {
		return
				dimension == DIMENSION::D3 &&
				info::ecf->structural_mechanics_3d.discretization.size() &&
				info::ecf->structural_mechanics_3d.discretization.begin()->second == PhysicsConfiguration::DISCRETIZATION::FEM_TDNNS;
	};

	LinearSystem *current = NULL;
	int DOFs = 0;
	Kernel *kernel = NULL, *pKernel = previous ? previous->assembler()->composer->kernel : NULL;
	switch (dimension) {
	case DIMENSION::D2: DOFs = 2;
		if (dynamic_cast<HeatTransfer2DKernel*>(pKernel)) {
			kernel = new StructuralMechanics2DKernel(dynamic_cast<HeatTransfer2DKernel*>(pKernel), physics, gsettings, loadStep); break;
		} else {
			kernel = new StructuralMechanics2DKernel(dynamic_cast<StructuralMechanics2DKernel*>(pKernel), physics, gsettings, loadStep); break;
		}

		break;
	case DIMENSION::D3: {
		DOFs = 3;
		if (tdnns()) {
			kernel = new StructuralMechanics3DTDNNSKernel(dynamic_cast<StructuralMechanics3DTDNNSKernel*>(pKernel), physics, gsettings, loadStep); break;
		} else {
			switch (loadStep.type) {
			case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
			case LoadStepSolverConfiguration::TYPE::TRANSIENT:
				kernel = new StructuralMechanics3DKernel(dynamic_cast<StructuralMechanics3DKernel*>(pKernel), physics, gsettings, loadStep); break;
			case LoadStepSolverConfiguration::TYPE::HARMONIC:
				if (dynamic_cast<HarmonicBalance3DKernel*>(pKernel)) {
					kernel = new HarmonicBalance3DKernel(dynamic_cast<HarmonicBalance3DKernel*>(pKernel), physics, gsettings, loadStep); break;
				} else {
					kernel = new HarmonicBalance3DKernel(dynamic_cast<StructuralMechanics3DKernel*>(pKernel), physics, gsettings, loadStep); break;
				}
			}
		}
	} break;
	default: break;
	}

	switch (loadStep.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI: {
		FETISystem *system = new FETISystem(1, 1, loadStep.feti);
		system->assembler()->composer = new NodesUniformFETIComposer(loadStep.feti, kernel, system->assembler(), DOFs);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE: {
		HYPRESystem *system = new HYPRESystem(1, 1, loadStep.hypre);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: {
		MKLPDSSSystem *system = new MKLPDSSSystem(1, 1, loadStep.mklpdss);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: {
		PARDISOSystem *system = new PARDISOSystem(1, 1, loadStep.pardiso);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: {
		SuperLUSystem *system = new SuperLUSystem(1, 1, loadStep.superlu);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::WSMP: {
		WSMPSystem *system = new WSMPSystem(1, 1, loadStep.wsmp);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, system->assembler(), DOFs);
		}
		current = system;
	} break;
	default:
		eslog::globalerror("Not implemented / unknown linear solver.\n");
	}

	switch (loadStep.type) {
	case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
	case LoadStepSolverConfiguration::TYPE::TRANSIENT:
		current->builder = new TimeBuilder(); break;
	case LoadStepSolverConfiguration::TYPE::HARMONIC:
		switch (dimension) {
		case DIMENSION::D2: current->builder = new HarmonicBuilder(2); break;
		case DIMENSION::D3: current->builder = new HarmonicBuilder(3); break;
		default: break;
		} break;
	default:
		eslog::globalerror("ESPRESO internal error: not implemented loadatep type.\n");
	}

	current->init();
	return current;
}

static SubStepSolver* getSubStepSolver(SubStepSolver *previous, HeatTransferLoadStepSolverConfiguration &loadStep, LinearSystem *system)
{
	SubStepSolver* current = NULL;
	switch (loadStep.mode) {
	case LoadStepSolverConfiguration::MODE::LINEAR: current = new LinearSubStep(system); break;
	case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new NewtonRaphson(system, loadStep.nonlinear_solver); break;
	default:
		eslog::globalerror("Not implemented sub-step solver.\n");
	}
	current->init(previous);
	return current;
}

static SubStepSolver* getSubStepSolver(SubStepSolver *previous, StructuralMechanicsLoadStepSolverConfiguration &loadStep, LinearSystem *system)
{
	SubStepSolver* current = NULL;
	switch (loadStep.mode) {
	case LoadStepSolverConfiguration::MODE::LINEAR: current = new LinearSubStep(system); break;
	case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new NewtonRaphson(system, loadStep.nonlinear_solver); break;
	default:
		eslog::globalerror("Not implemented sub-step solver.\n");
	}
	current->init(previous);
	return current;
}

static LoadStepSolver* getLoadStepSolver(LoadStepSolver *previous, HeatTransferLoadStepConfiguration &loadStep, LinearSystem *system, SubStepSolver *subStepSolver)
{
	LoadStepSolver* current = NULL;
	if (loadStep.topology_optimization) {
		current = new TopologyOptimization(system, subStepSolver, loadStep.topology_optimization_settings);
	} else {
		switch (loadStep.type) {
		case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
			switch (loadStep.mode){
			case LoadStepSolverConfiguration::MODE::LINEAR: current = new SteadyStateSolver(system, subStepSolver, loadStep.duration_time); break;
			case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new PseudoTimeStepping(system, subStepSolver, loadStep.nonlinear_solver, loadStep.duration_time); break;
			} break;
		case LoadStepSolverConfiguration::TYPE::TRANSIENT: current = new TransientFirstOrderImplicit(system, subStepSolver, loadStep.transient_solver, loadStep.duration_time); break;
		default:
			eslog::globalerror("Not implemented load-step solver.\n");
		}
	}

	current->init(previous);
	return current;
}

static LoadStepSolver* getLoadStepSolver(LoadStepSolver *previous, StructuralMechanicsLoadStepConfiguration &loadStep, LinearSystem *system, SubStepSolver *subStepSolver)
{
	LoadStepSolver* current = NULL;
	if (loadStep.topology_optimization) {
		current = new TopologyOptimization(system, subStepSolver, loadStep.topology_optimization_settings);
	} else {
		switch (loadStep.type) {
		case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
			switch (loadStep.mode){
			case LoadStepSolverConfiguration::MODE::LINEAR: current = new SteadyStateSolver(system, subStepSolver, loadStep.duration_time); break;
			case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new PseudoTimeStepping(system, subStepSolver, loadStep.nonlinear_solver, loadStep.duration_time); break;
			} break;
		case LoadStepSolverConfiguration::TYPE::TRANSIENT: current = new TransientSecondOrderImplicit(system, subStepSolver, loadStep.transient_solver, loadStep.duration_time); break;
		case LoadStepSolverConfiguration::TYPE::HARMONIC: current = new HarmonicSolver(system, subStepSolver, loadStep.harmonic_solver, loadStep.duration_time); break;
		default:
			eslog::globalerror("Not implemented load-step solver.\n");
		}
	}

	current->init(previous);
	return current;
}

void PhysicalSolver::run()
{
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D:
		{ PhysicalSolver first; PhysicalSolver second; runCoupled(first, second, info::ecf->thermo_elasticity_2d); } break;
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D:
		{ PhysicalSolver first; PhysicalSolver second; runCoupled(first, second, info::ecf->thermo_elasticity_3d); } break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		{ PhysicalSolver solver; runSingle<>(solver, info::ecf->heat_transfer_2d); } break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		{ PhysicalSolver solver; runSingle(solver, info::ecf->heat_transfer_3d); } break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		{ PhysicalSolver solver; runSingle(solver, info::ecf->structural_mechanics_2d); } break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		{ PhysicalSolver solver; runSingle(solver, info::ecf->structural_mechanics_3d); } break;
	default:
		eslog::globalerror("Physical solver: not implemented physical solver.\n");
	}
}

PhysicalSolver::PhysicalSolver()
: loadStepSolver(NULL), subStepSolver(NULL), system(NULL)
{

}

PhysicalSolver::~PhysicalSolver()
{
	clear();
}

void PhysicalSolver::clear()
{
	if (loadStepSolver) { delete loadStepSolver; }
	if (subStepSolver) { delete subStepSolver; }
	if (system) { delete system; }
}

template <typename TPhysics>
void PhysicalSolver::runSingle(PhysicalSolver &solver, TPhysics &configuration)
{
	while (step::loadstep < configuration.load_steps) {
		eslog::nextStep(step::loadstep + 1);

		eslog::startln("PHYSICS BUILDER: LOAD STEP STARTED", "PHYSICS BUILDER");

		auto &loadStepSettings = configuration.load_steps_settings.at(step::loadstep + 1);

		PhysicalSolver prev = solver;
		solver.system = getSystem(prev.system, configuration, configuration, loadStepSettings, configuration.dimension);
		eslog::checkpointln("PHYSICS BUILDER: LINEAR SYSTEM COMPOSED");

		solver.subStepSolver = getSubStepSolver(prev.subStepSolver, loadStepSettings, solver.system);
		eslog::checkpointln("PHYSICS BUILDER: SUBSTEP SOLVER INITTED");

		solver.loadStepSolver = getLoadStepSolver(prev.loadStepSolver, loadStepSettings, solver.system, solver.subStepSolver);
		eslog::checkpointln("PHYSICS BUILDER: LOAD STEP SOLVER INITTED");

		switch(info::ecf->mesh_morphing.type) {
		case MORPHING_TYPE::NONE: break;
		case MORPHING_TYPE::RBF: morphing::rbf(); break;
		}

		if (step::isInitial()) {
			info::mesh->store->updateMonitors();
		}

		eslog::endln("PHYSICS BUILDER: MONITORS SET");

		eslog::checkpoint("ESPRESO: PHYSICS PREPARED");
		eslog::param("LOADSTEP", step::loadstep + 1);
		eslog::ln();

		eslog::startln("PHYSICS SOLVER: STARTED", "PHYSICS SOLVER");
		solver.loadStepSolver->run();
		eslog::endln("PHYSICS SOLVER: FINISHED");

		++step::loadstep;
		eslog::checkpoint("ESPRESO: PHYSICS SOLVED");
		eslog::param("LOADSTEP", step::loadstep);
		eslog::ln();
	}
}

template <typename TPhysics>
void PhysicalSolver::runCoupled(PhysicalSolver &first, PhysicalSolver &second, TPhysics &configuration)
{
	while (step::loadstep < configuration.load_steps) {
		eslog::nextStep(step::loadstep + 1);

		eslog::startln("PHYSICS BUILDER: LOAD STEP STARTED", "PHYSICS BUILDER");
		auto &loadStepSettings = configuration.load_steps_settings.at(step::loadstep + 1);

		PhysicalSolver _first = first, _second = second;
		first.system = getSystem(_first.system, configuration, configuration, loadStepSettings.heat_transfer, configuration.dimension);
		second.system = getSystem(first.system, configuration, configuration, loadStepSettings.structural_mechanics, configuration.dimension);
		eslog::checkpointln("PHYSICS BUILDER: LINEAR SYSTEM COMPOSED");

		first.subStepSolver = getSubStepSolver(_first.subStepSolver, loadStepSettings.heat_transfer, first.system);
		second.subStepSolver = getSubStepSolver(_second.subStepSolver, loadStepSettings.structural_mechanics, second.system);
		eslog::checkpointln("PHYSICS BUILDER: SUBSTEP SOLVER INITTED");

		first.loadStepSolver = getLoadStepSolver(_first.loadStepSolver, loadStepSettings.heat_transfer, first.system, first.subStepSolver);
		second.loadStepSolver = getLoadStepSolver(_second.loadStepSolver, loadStepSettings.structural_mechanics, second.system, second.subStepSolver);
		eslog::checkpointln("PHYSICS BUILDER: LOAD STEP SOLVER INITTED");

		if (step::isInitial()) {
			info::mesh->store->updateMonitors();
		}
		eslog::endln("PHYSICS BUILDER: MONITORS SET");

		eslog::checkpoint("ESPRESO: PHYSICS PREPARED");
		eslog::param("LOADSTEP", step::loadstep + 1);
		eslog::ln();

		eslog::startln("PHYSICS SOLVER: STARTED", "PHYSICS SOLVER");
		step::substep = step::duplicate::offset;
		step::iteration = 0;

		while (!step::isLast()) {
			info::mesh->store->suppress();
			first.loadStepSolver->runNextSubstep();
			info::mesh->store->permit();
			second.loadStepSolver->runNextSubstep();
			step::substep++;
			step::iteration = 0;
		}
		eslog::endln("PHYSICS SOLVER: FINISHED");

		++step::loadstep;
		eslog::checkpoint("ESPRESO: PHYSICS SOLVED");
		eslog::param("LOADSTEP", step::loadstep);
		eslog::ln();
	}
}


