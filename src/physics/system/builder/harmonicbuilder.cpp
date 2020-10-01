
#include "harmonicbuilder.h"
#include "esinfo/stepinfo.h"
#include "physics/system/wsmpsystem.h"
#include "physics/system/superlusystem.h"
#include "physics/system/pardisosystem.h"
#include "physics/system/mklpdsssystem.h"
#include "physics/system/hypresystem.h"
#include "physics/system/fetisystem.h"
#include "physics/composer/composer.h"
#include "physics/composer/feti/feti.composer.h"
#include "physics/kernels/kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "config/ecf/linearsolver/feti.h"

#include <cmath>

using namespace espreso;

void HarmonicBuilder::init(AssemblerData &assembler, SolverData &solver)
{
	solver.K->uniformCombination(assembler.K, assembler.K, DOFs, DOFs);
	solver.f->uniformCombination(assembler.f->at(0), assembler.f->at(1), DOFs, DOFs);
	solver.R->shallowCopyStructure(solver.f);
	solver.x->shallowCopyStructure(solver.f);
	solver.y->shallowCopyStructure(solver.f);
	solver.BC->uniformCombination(assembler.BC->at(0), assembler.BC->at(0), DOFs, DOFs); // BC does not have 2 vectors

	solver.x->at(0)->fillCombinedValues(&assembler.composer->kernel->solutions[0], 0, DOFs, 2 * DOFs);
	solver.x->at(0)->fillCombinedValues(&assembler.composer->kernel->solutions[1], DOFs, DOFs, 2 * DOFs);

	reset(Builder::Request::KCM | Builder::Request::RBCf, assembler, solver);
}

void HarmonicBuilder::buildSystem(AssemblerData &assembler, SolverData &solver)
{
	if (rayleighDamping || coriolisDamping) {
		solver.K->type = MatrixType::REAL_UNSYMMETRIC;
	} else {
		solver.K->type = MatrixType::REAL_SYMMETRIC_INDEFINITE;
	}

	if (step::type == step::TYPE::FREQUENCY) {
		if (matrices & Builder::Request::K) {
			solver.K->addToCombination(1, assembler.K, 0, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(1, assembler.K, DOFs, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
		}
		if (matrices & Builder::Request::C) {
			solver.K->addToCombination(-step::frequency::angular, assembler.C, 0, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination( step::frequency::angular, assembler.C, DOFs, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(1, assembler.CM, 0, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(1, assembler.CM, DOFs, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
		} else {
			if (rayleighDamping) {
				double stiffCoef = stiffnessDamping + structuralDampingCoefficient / step::frequency::angular;
				solver.K->addToCombination(-step::frequency::angular * stiffCoef, assembler.K, 0, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
				solver.K->addToCombination(-step::frequency::angular * massDamping, assembler.M, 0, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
				solver.K->addToCombination( step::frequency::angular * stiffCoef, assembler.K, DOFs, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
				solver.K->addToCombination( step::frequency::angular * massDamping, assembler.M, DOFs, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			}
		}
		if (matrices & Builder::Request::M) {
			solver.K->addToCombination(-step::frequency::angular * step::frequency::angular, assembler.M, 0, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(-step::frequency::angular * step::frequency::angular, assembler.M, DOFs, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
		}

		if (matrices & Builder::Request::R) {
			solver.R->at(0)->addToCombination(1, assembler.R->at(0), 0, DOFs, 2 * DOFs);
			solver.R->at(0)->addToCombination(1, assembler.R->at(1), DOFs, DOFs, 2 * DOFs);
		}

		solver.f->at(0)->fillCombinedValues(assembler.f->at(0), 0, DOFs, 2 * DOFs);
		solver.f->at(0)->fillCombinedValues(assembler.f->at(1), DOFs, DOFs, 2 * DOFs);
		solver.BC->at(0)->fillCombinedValues(assembler.BC->at(0), 0, DOFs, 2 * DOFs);
		solver.BC->at(0)->fillCombinedValues(assembler.BC->at(0), DOFs, DOFs, 2 * DOFs);
	}

	if (step::type == step::TYPE::FTT) {
//		if (matrices & Builder::Request::M) {
//			solver.K->addToCombination(-step::frequency::angular * step::frequency::angular, assembler.M, 0, 0, DOFs, 2 * DOFs);
//			solver.K->addToCombination(-step::frequency::angular * step::frequency::angular, assembler.M, DOFs, DOFs, DOFs, 2 * DOFs);
//		}
//		if (matrices & Builder::Request::C) {
//			solver.K->addToCombination(-1, assembler.C, 0, DOFs, DOFs, 2 * DOFs);
//			solver.K->addToCombination( 1, assembler.C, DOFs, 0, DOFs, 2 * DOFs);
//		}

		double cosVal = std::cos(step::ftt::step * 2 * M_PI);
		double sinVal = std::sin(step::ftt::step * 2 * M_PI);

		double coeff = 2.0 / step::ftt::steps;

		if (matrices & Builder::Request::K) {
			solver.K->addToCombination(cosVal * cosVal * coeff, assembler.K, 0, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(sinVal * sinVal * coeff, assembler.K, DOFs, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(cosVal * sinVal * coeff, assembler.K, 0, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			solver.K->addToCombination(sinVal * cosVal * coeff, assembler.K, DOFs, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
		}
		if (matrices & Builder::Request::R) {
			solver.R->at(0)->addToCombination(cosVal * coeff, assembler.R->at(0), 0, DOFs, 2 * DOFs);
			solver.R->at(0)->addToCombination(sinVal * coeff, assembler.R->at(1), DOFs, DOFs, 2 * DOFs);
		}
	}
}

void HarmonicBuilder::updateSolution(AssemblerData &assembler, SolverData &solver)
{
	assembler.x->at(0)->fillValuesFromCombination(solver.x->at(0), 0, DOFs, 2 * DOFs);
	assembler.x->at(1)->fillValuesFromCombination(solver.x->at(0), DOFs, DOFs, 2 * DOFs);
	for (esint n = 0; n < assembler.x->nvectors; n++) {
		assembler.composer->kernel->solutions[n].fillData(assembler.x->at(n));
	}
	assembler.composer->kernel->solutionChanged();
}

void HarmonicBuilder::init(WSMPSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::init(SuperLUSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::init(PARDISOSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::init(MKLPDSSSystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::init(HYPRESystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
	system.assemblers[0].numFnc = system.assemblers[0].composer->kernel->solverDataProvider->hypre->numfnc();
//	system.assemblers[0].composer->assembler->initHypreKernels(system.assemblers[0].K, system.assemblers[0].N);
//	system.solvers[0].N.uniformCombination(&system.assemblers[0].N, &system.assemblers[0].N, DOFs, DOFs);
}

void HarmonicBuilder::init(FETISystem &system)
{
	init(system.assemblers[0], system.solvers[0]);
	system.solvers[0].mortars.shallowCopy(&system.assemblers[0].mortars);
	system.solvers[0].buildB1();
	system.solvers[0].Kdiag.shallowCopyStructure(system.solvers[0].f.at(0));
	if (system.solvers[0].solver.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
		system.assemblers[0].composer->kernel->solverDataProvider->feti->initKernels(
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
			system.solvers[0].B0.uniformCombination(&system.assemblers[0].B0, &system.assemblers[0].B0, DOFs, DOFs);
			system.solvers[0].B0.addToCombination(1, &system.assemblers[0].B0, 0, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
			system.solvers[0].B0.addToCombination(1, &system.assemblers[0].B0, DOFs, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
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

void HarmonicBuilder::buildSystem(WSMPSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::buildSystem(SuperLUSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::buildSystem(PARDISOSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::buildSystem(MKLPDSSSystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::buildSystem(HYPRESystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
	system.solvers[0].numFnc = 2 * system.assemblers[0].numFnc;
//	system.assemblers[0].composer->assembler->fillHypreKernels(system.assemblers[0].K, system.assemblers[0].N);
//	system.solvers[0].N.fillCombinedValues(&system.assemblers[0].N, 0, DOFs, 2 * DOFs);
//	system.solvers[0].N.fillCombinedValues(&system.assemblers[0].N, DOFs, DOFs, 2 * DOFs);
}

void HarmonicBuilder::buildSystem(FETISystem &system)
{
	buildSystem(system.assemblers[0], system.solvers[0]);
	if (system.solvers[0].solver.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
		system.assemblers[0].composer->kernel->solverDataProvider->feti->fillKernels(
				system.assemblers[0].K, system.assemblers[0].M, system.assemblers[0].N1, system.assemblers[0].N2, system.assemblers[0].RegMat,
				system.solvers[0].solver.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
		system.solvers[0].N1.fillData(&system.assemblers[0].N1);
		system.solvers[0].N2.fillData(&system.assemblers[0].N2);
		system.solvers[0].RegMat.fillData(&system.assemblers[0].RegMat);

//		int NDIM = system.solvers[0].kernelDimension;
//		system.solvers[0].N1.fill(0);
//		system.solvers[0].N1.addToCombination(1, &system.assemblers[0].N1, 0, 0, DOFs, NDIM, 2 * DOFs, 2 * NDIM);
//		system.solvers[0].N1.addToCombination(1, &system.assemblers[0].N1, DOFs, NDIM, DOFs, NDIM, 2 * DOFs, 2 * NDIM);
//		system.solvers[0].N2.fill(0);
//		system.solvers[0].N2.addToCombination(1, &system.assemblers[0].N2, 0, 0, DOFs, NDIM, 2 * DOFs, 2 * NDIM);
//		system.solvers[0].N2.addToCombination(1, &system.assemblers[0].N2, DOFs, NDIM, DOFs, NDIM, 2 * DOFs, 2 * NDIM);
//		system.solvers[0].RegMat.fill(0);
//		system.solvers[0].RegMat.addToCombination(1, &system.assemblers[0].RegMat, 0, 0, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
//		system.solvers[0].RegMat.addToCombination(1, &system.assemblers[0].RegMat, DOFs, DOFs, DOFs, DOFs, 2 * DOFs, 2 * DOFs);
	}
	for (size_t i = 0; i < system.solvers.size(); i++) {
		if (
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
				system.solvers[i].solver.configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

			system.solvers[i].origK.fillData(&system.solvers[i].K);
		}
	}
}

void HarmonicBuilder::updateSolution(WSMPSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::updateSolution(SuperLUSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::updateSolution(PARDISOSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::updateSolution(MKLPDSSSystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::updateSolution(HYPRESystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

void HarmonicBuilder::updateSolution(FETISystem &system)
{
	updateSolution(system.assemblers[0], system.solvers[0]);
}

