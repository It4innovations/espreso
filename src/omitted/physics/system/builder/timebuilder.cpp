
#include "timebuilder.h"

#include "physics/system/wsmpsystem.h"
#include "physics/system/superlusystem.h"
#include "physics/system/pardisosystem.h"
#include "physics/system/mklpdsssystem.h"
#include "physics/system/hypresystem.h"
#include "physics/system/fetisystem.h"
#include "physics/composer/composer.h"
#include "physics/composer/feti/feti.composer.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "config/ecf/linearsolver/feti.h"

using namespace espreso;

void TimeBuilder::init(AssemblerData &assembler, SolverData &solver)
{
	solver.K->shallowCopyStructure(assembler.K);
	solver.f->shallowCopyStructure(assembler.f);
	solver.R->shallowCopyStructure(assembler.R);
	solver.x->shallowCopyStructure(assembler.x);
	solver.y->shallowCopyStructure(assembler.x);
	solver.BC->shallowCopyStructure(assembler.BC);

	solver.x->at(0)->fillData(assembler.x->at(0));
}

void TimeBuilder::buildSystem(AssemblerData &assembler, SolverData &solver)
{
	if (matrices & Builder::Request::K) {
		solver.K->add(timeIntegrationConstantK, assembler.K);
	}
	if (matrices & Builder::Request::C) {
		solver.K->add(timeIntegrationConstantC, assembler.C);
	}
	if (matrices & Builder::Request::M) {
		solver.K->add(timeIntegrationConstantM, assembler.M);
	}

	if (matrices & Builder::Request::f) {
		solver.f->fillData(assembler.f);
	}
	if (matrices & Builder::Request::R) {
		solver.R->fillData(assembler.R);
	}
	if (matrices & Builder::Request::BC) {
		solver.BC->fillData(assembler.BC);
	}
}

void TimeBuilder::updateSolution(AssemblerData &assembler, SolverData &solver)
{
	assembler.x->fillData(solver.x);
	assembler.composer->solutionChanged(assembler.x);
}

void TimeBuilder::init(WSMPSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::init(SuperLUSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::init(PARDISOSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::init(MKLPDSSSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::init(HYPRESystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
	system.assemblers[0].numFnc = system.assemblers[0].composer->provider()->hypre->numfnc();
	system.assemblers[0].composer->provider()->hypre->initKernels(system.assemblers[0].K, system.assemblers[0].N);
	system.solvers[0].N.shallowCopyStructure(&system.assemblers[0].N);
}

void TimeBuilder::init(FETISystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
	system.solvers[0].mortars.shallowCopy(&system.assemblers[0].mortars);
	system.solvers[0].gapDirection.shallowCopy(&system.assemblers[0].gapDirection);
	system.solvers[0].gap.shallowCopy(&system.assemblers[0].gap);
	system.solvers[0].buildB1();
	system.solvers[0].Kdiag.shallowCopyStructure(system.solvers[0].f.at(0));
	if (system.solvers[0].solver.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
		system.assemblers[0].composer->provider()->feti->initKernels(
				system.assemblers[0].K, system.assemblers[0].N1, system.assemblers[0].N2, system.assemblers[0].RegMat,
				system.solvers[0].solver.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
		system.solvers[0].N1.shallowCopyStructure(&system.assemblers[0].N1);
		system.solvers[0].N2.shallowCopyStructure(&system.assemblers[0].N2);
		system.solvers[0].RegMat.shallowCopyStructure(&system.assemblers[0].RegMat);
	}
	if (system.solvers[0].solver.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI) {
		switch (system.solvers[0].solver.configuration.B0_type) {
		case FETIConfiguration::B0_TYPE::CORNERS:
			reinterpret_cast<FETIComposer*>(system.assemblers[0].composer)->buildB0FromCorners(system.assemblers[0].B0);
			system.solvers[0].B0.shallowCopy(&system.assemblers[0].B0);
			break;
		case FETIConfiguration::B0_TYPE::KERNELS:
			system.solvers[0].buildB0();
			break;
		}
	}
	for (size_t i = 0; i < system.solvers.size(); i++) {
		if (
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

			system.solvers[i].origK.shallowCopyStructure(&system.solvers[i].K);
		}
	}
}

void TimeBuilder::buildSystem(WSMPSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::buildSystem(SuperLUSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::buildSystem(PARDISOSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::buildSystem(MKLPDSSSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::buildSystem(HYPRESystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
	system.solvers[0].numFnc = system.assemblers[0].numFnc;
	system.assemblers[0].composer->provider()->hypre->fillKernels(system.assemblers[0].K, system.assemblers[0].N);
	system.solvers[0].N.fillData(&system.assemblers[0].N);
}

void TimeBuilder::buildSystem(FETISystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
	if (system.solvers[0].solver.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
		system.assemblers[0].composer->provider()->feti->fillKernels(
				system.assemblers[0].K, system.assemblers[0].M, system.assemblers[0].N1, system.assemblers[0].N2, system.assemblers[0].RegMat,
				system.solvers[0].solver.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
		system.solvers[0].N1.fillData(&system.assemblers[0].N1);
		system.solvers[0].N2.fillData(&system.assemblers[0].N2);
		system.solvers[0].RegMat.fillData(&system.assemblers[0].RegMat);
	}
	for (size_t i = 0; i < system.solvers.size(); i++) {
		if (
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

			system.solvers[i].origK.fillData(&system.solvers[i].K);
		}
	}
}

void TimeBuilder::updateSolution(WSMPSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::updateSolution(SuperLUSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::updateSolution(PARDISOSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::updateSolution(MKLPDSSSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::updateSolution(HYPRESystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void TimeBuilder::updateSolution(FETISystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

