
#include "loadstepiterator.h"

#include "basis/expression/expression.h"
#include "basis/evaluator/expressionevaluator.h"
#include "basis/utilities/utils.h"
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
#include "physics/assembler/modules/heattransfer.module.opt.h"
#include "physics/kernels/heattransfer2d.kernel.h"
#include "physics/kernels/heattransfer3d.kernel.h"
#include "physics/kernels/structuralmechanics2d.kernel.h"
#include "physics/kernels/structuralmechanics3d.elasticity.kernel.h"
#include "physics/kernels/structuralmechanics3d.harmonic.kernel.h"
#include "physics/kernels/structuralmechanics3d.tdnns.kernel.h"
#include "output/output.h"

#include <algorithm>

using namespace espreso;

static LinearSystem* getSystem(LinearSystem *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	LinearSystem *current = NULL;
	Kernel *kernel = NULL, *pKernel = previous ? previous->assembler()->composer->kernel : NULL;
	ModuleOpt *opt = NULL;
	switch (gsettings.kernel) {
	case HeatTransferGlobalSettings::KERNEL::OLD:
		switch (dimension) {
		case DIMENSION::D2: kernel = new HeatTransfer2DKernel(dynamic_cast<HeatTransfer2DKernel*>(pKernel), physics, gsettings, loadStep); break;
		case DIMENSION::D3: kernel = new HeatTransfer3DKernel(dynamic_cast<HeatTransfer3DKernel*>(pKernel), physics, gsettings, loadStep); break;
		default: break;
		}
		break;
	case HeatTransferGlobalSettings::KERNEL::OPT:
	case HeatTransferGlobalSettings::KERNEL::VEC:
		opt = new HeatTransferModuleOpt(dynamic_cast<HeatTransferModuleOpt*>(pKernel), gsettings, loadStep);
		break;
	}

	switch (loadStep.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI: {
		FETISystem *system = new FETISystem(1, 1, loadStep.feti);
		system->assembler()->composer = new NodesUniformFETIComposer(loadStep.feti, kernel, opt, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE: {
		HYPRESystem *system = new HYPRESystem(1, 1, loadStep.hypre);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, opt, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: {
		MKLPDSSSystem *system = new MKLPDSSSystem(1, 1, loadStep.mklpdss);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, opt, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: {
		PARDISOSystem *system = new PARDISOSystem(1, 1, loadStep.pardiso);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, opt, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: {
		SuperLUSystem *system = new SuperLUSystem(1, 1, loadStep.superlu);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, opt, system->assembler(), 1);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::WSMP: {
		WSMPSystem *system = new WSMPSystem(1, 1, loadStep.wsmp);
		system->assembler()->composer = new NodesUniformDistributedComposer(kernel, opt, system->assembler(), 1);
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
	case DIMENSION::D2: DOFs = 2; kernel = new StructuralMechanics2DKernel(dynamic_cast<StructuralMechanics2DKernel*>(pKernel), physics, gsettings, loadStep); break;
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
		system->assembler()->composer = new NodesUniformFETIComposer(loadStep.feti, kernel, NULL, system->assembler(), DOFs);
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE: {
		HYPRESystem *system = new HYPRESystem(1, 1, loadStep.hypre);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, NULL, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, NULL, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: {
		MKLPDSSSystem *system = new MKLPDSSSystem(1, 1, loadStep.mklpdss);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, NULL, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, NULL, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: {
		PARDISOSystem *system = new PARDISOSystem(1, 1, loadStep.pardiso);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, NULL, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, NULL, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: {
		SuperLUSystem *system = new SuperLUSystem(1, 1, loadStep.superlu);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, NULL, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, NULL, system->assembler(), DOFs);
		}
		current = system;
	} break;
	case LoadStepSolverConfiguration::SOLVER::WSMP: {
		WSMPSystem *system = new WSMPSystem(1, 1, loadStep.wsmp);
		if (tdnns()) {
			system->assembler()->composer = new FacesEdgesUniformDistributedComposer(kernel, NULL, system->assembler(), 4, 2);
		} else {
			system->assembler()->composer = new NodesUniformDistributedComposer(kernel, NULL, system->assembler(), DOFs);
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
		case LoadStepSolverConfiguration::TYPE::HARMONIC: current = new HarmonicSolver(system, subStepSolver, loadStep, loadStep.duration_time); break;
		default:
			eslog::globalerror("Not implemented load-step solver.\n");
		}
	}

	current->init(previous);
	return current;
}


LoadStepIterator::LoadStepIterator()
: _loadStepSolver(NULL), _subStepSolver(NULL), _system(NULL)
{

}

LoadStepIterator::~LoadStepIterator()
{
	if (_loadStepSolver != NULL) { delete _loadStepSolver; }
	if (_subStepSolver != NULL) { delete _subStepSolver; }
	if (_system != NULL) { delete _system; }
}

void LoadStepIterator::prepareExpressions()
{
	if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D && info::ecf->heat_transfer_2d.kernel == HeatTransferConfiguration::KERNEL::OLD) {
		return;
	}
	if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D && info::ecf->heat_transfer_3d.kernel == HeatTransferConfiguration::KERNEL::OLD) {
		return;
	}
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D: HeatTransferModuleOpt::createParameters(); break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D: HeatTransferModuleOpt::createParameters(); break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
	default: eslog::globalerror("Unknown physics.\n");
	}

	struct __param__ {
		std::string value, replacement;
		std::vector<std::string> parameters;

		__param__(const std::string &param): value(param), replacement(value), parameters({param}) {}
		__param__(const std::string &value, const std::vector<std::string> &parameters)
		: value(value), parameters(parameters)
		{
			if (parameters.size()) {
				replacement = "sqrt(";
				for (auto p = parameters.begin(); p != parameters.end(); ++p) {
					if (p != parameters.begin()) {
						replacement += " + ";
					}
					replacement += *p + " * " + *p;
				}
				replacement += ")";
			} else {
				this->replacement = value;
				this->parameters.push_back(value);
			}
		}
	};

	std::vector<__param__> parameters = { {"TIME"} };
	if (info::mesh->dimension == 2) {
		parameters.push_back({ "COORDINATE_X" });
		parameters.push_back({ "COORDINATE_Y" });
		parameters.push_back({ "COORDINATE", { "COORDINATE_X", "COORDINATE_Y" }});
	}
	if (info::mesh->dimension == 3) {
		parameters.push_back({ "COORDINATE_X" });
		parameters.push_back({ "COORDINATE_Y" });
		parameters.push_back({ "COORDINATE_Z" });
		parameters.push_back({ "COORDINATE", { "COORDINATE_X", "COORDINATE_Y", "COORDINATE_Z" }});
	}

	std::vector<const NamedData*> data;
	data.insert(data.end(), info::mesh->nodes->data.begin(), info::mesh->nodes->data.end());
	data.insert(data.end(), info::mesh->elements->data.begin(), info::mesh->elements->data.end());
	for (auto it = data.begin(); it != data.end(); ++it) {
		const NamedData *p = *it;
		if (p->dataType == NamedData::DataType::SCALAR) {
			parameters.push_back({p->name});
		} else {
			std::vector<std::string> params;
			for (int d = 0; d < p->dimension; ++d) {
				params.push_back(p->name + p->suffix(d));
				parameters.push_back(params.back());
			}
			parameters.push_back({ p->name, params });
		}
	}

//	for (auto fnc = info::ecf->functions.begin(); fnc != info::ecf->functions.end(); ++fnc) {
//		parameters.push_back({ fnc->first });
//	}

	for (auto it = ECFExpression::parametrized.begin(); it != ECFExpression::parametrized.end(); ++it) {
		ECFExpression *expression = *it;
		std::string upper = expression->value;
		for (size_t i = 0; i < upper.size(); i++) {
			upper[i] = std::toupper(upper[i]);
		}

		size_t pindex = 0;
		std::vector<std::pair<size_t, size_t> > pindices;
		for (auto param = parameters.begin(); param != parameters.end(); ++param, ++pindex) {
			size_t i = 0;
			while ((i = upper.find(param->value, i)) != std::string::npos) {
				pindices.push_back(std::make_pair(i++, pindex));
			}
		}
		std::sort(pindices.begin(), pindices.end());

		for (size_t i = 1; i < pindices.size(); ++i) {
			if (pindices[i - 1].first == pindices[i].first) {
				if (parameters[pindices[i].second].value.size() < parameters[pindices[i - 1].second].value.size()) {
					pindices[i] = pindices[i - 1];
				}
				pindices[i - 1].first = std::string::npos;
			}
		}

		for (size_t i = 0, offset = 0; i < pindices.size(); ++i) {
			if (pindices[i].first != std::string::npos) {
				expression->parameters.insert(expression->parameters.end(), parameters[pindices[i].second].parameters.begin(), parameters[pindices[i].second].parameters.end());
				expression->value.replace(pindices[i].first + offset, parameters[pindices[i].second].value.size(), parameters[pindices[i].second].replacement);
				offset += parameters[pindices[i].second].replacement.size() - parameters[pindices[i].second].value.size();
			}
		}

		utils::sortAndRemoveDuplicates(expression->parameters);

		if (Expression::isValid(expression->value, expression->parameters)) {
			expression->evaluator = new ExpressionEvaluator(expression->value, expression->parameters);
		} else {
			eslog::globalerror("PARSE ERROR: Expression '%s' cannot be evaluated.\n", expression->value.c_str());
		}
	}

	// TODO: check circular dependency
}

bool LoadStepIterator::next()
{
	eslog::nextStep(step::step.loadstep + 1);
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		return next(info::ecf->heat_transfer_2d);
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		return next(info::ecf->heat_transfer_3d);
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		return next(info::ecf->structural_mechanics_2d);
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		return next(info::ecf->structural_mechanics_3d);
	default:
		eslog::globalerror("Unknown physics.\n");
	}

	return false;
}

template <typename TPhysics>
bool LoadStepIterator::next(TPhysics &configuration)
{
	eslog::startln("PHYSICS BUILDER: LOAD STEP STARTED", "PHYSICS BUILDER");

	auto &loadStepSettings = configuration.load_steps_settings.at(step::step.loadstep + 1);
	if (loadStepSettings.topology_optimization && loadStepSettings.type != LoadStepSolverConfiguration::TYPE::STEADY_STATE) {
		eslog::error("ESPRESO: not implemented feature: topology optimization can be used only for the STEADY STATE solver type.\n");
	}

	auto system = getSystem(_system, configuration, configuration, loadStepSettings, configuration.dimension);
	eslog::checkpointln("PHYSICS BUILDER: LINEAR SYSTEM COMPOSED");

	auto subStepSolver = getSubStepSolver(_subStepSolver, loadStepSettings, system);
	eslog::checkpointln("PHYSICS BUILDER: SUBSTEP SOLVER INITTED");

	auto loadStepSolver = getLoadStepSolver(_loadStepSolver, loadStepSettings, system, subStepSolver);
	eslog::checkpointln("PHYSICS BUILDER: LOAD STEP SOLVER INITTED");

	if (_system) { delete _system; }
	if (_subStepSolver) { delete _subStepSolver; }
	if (_loadStepSolver) { delete _loadStepSolver; }

	_system = system;
	_subStepSolver = subStepSolver;
	_loadStepSolver = loadStepSolver;

	if (step::isInitial()) {
		info::mesh->output->updateMonitors(step::step.type);
	}
	eslog::endln("PHYSICS BUILDER: MONITORS SET");

	eslog::checkpoint("ESPRESO: PHYSICS PREPARED");
	eslog::param("LOADSTEP", step::step.loadstep + 1);
	eslog::ln();

	eslog::startln("PHYSICS SOLVER: STARTED", "PHYSICS SOLVER");
	_loadStepSolver->run();
	eslog::endln("PHYSICS SOLVER: FINISHED");

	return ++step::step.loadstep < configuration.load_steps;
}
