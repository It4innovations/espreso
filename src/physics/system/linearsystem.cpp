
#include "linearsystem.h"
#include "builder/builder.h"
#include "output/output.h"
#include "physics/composer/composer.h"
#include "physics/kernels/kernel.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"
#include "math/matrix.h"
#include "math/vector.h"
#include "esinfo/autoopt.h"

#include <cstddef>

using namespace espreso;

AssemblerData::AssemblerData(Matrix *K, Matrix *M, Matrix *C, Matrix *CM, Vectors *R, Vectors *f, Vectors *x, Vectors *BC)
: K(K), M(M), C(C), CM(CM), R(R), f(f), x(x), BC(BC),
  composer(NULL)
{

}

AssemblerData::~AssemblerData()
{
	if (composer != NULL) { delete composer; }
}

SolverData::SolverData(Matrix *K, Vectors *R, Vectors *f, Vectors *x, Vectors *y, Vectors *BC, SystemSolver *linearSolver)
: K(K), R(R), f(f), x(x), y(y), BC(BC),
  linearSolver(linearSolver)
{

}

SolverData::~SolverData()
{
//	if (provider != NULL) { delete provider; }
//	if (linearSolver != NULL) { delete linearSolver; }
}

LinearSystem::~LinearSystem()
{
	if (builder) { delete builder; }
}

void LinearSystem::init()
{
	assembler()->composer->init();
	_builderInit();
	for (esint i = 0; i < nsolvers(); ++i) {
		solver(i)->linearSolver->init();
	}
}

void LinearSystem::nextSubstep()
{
	assembler()->composer->kernel->nextSubstep();
	eslog::checkpointln("PHYSICS SOLVER: PARAMETERS EVALUATED");
}

void LinearSystem::assemble()
{
	_builderReset();

	assembler()->composer->assemble(*builder);
	_builderCreateSystem();

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE ASSEMBLED MATRICES\n");
		if (nassemblers() == 1) {
			std::string prefix = utils::debugDirectory() + "/assembler";
			assembler()->print(builder, prefix.c_str(), "");
		} else {
			for (int i = 0; i < nassemblers(); i++) {
				std::string prefix = utils::debugDirectory() + "/assembler" + std::to_string(i + 1);
				assembler()->print(builder, prefix.c_str(), "");
			}
		}
	}

	if (step::step.type == step::TYPE::FREQUENCY && builder->AFTSamples && (builder->matrices & Builder::Request::R)) {
		step::step.type = step::TYPE::FTT;
		step::ftt.steps = builder->AFTSamples;
		step::ftt.period = 1 / step::frequency.current;
		for (step::ftt.step = 0; step::ftt.step < step::ftt.steps; ++step::ftt.step) {
			step::ftt.time = (double)step::ftt.step / step::ftt.steps;

			assembler()->composer->kernel->processSolution();
			assembler()->composer->kernel->solutionChanged();
			assembler()->composer->assemble(*builder);
			_builderCreateSystem();

			if (info::ecf->output.print_matrices) {
				eslog::storedata(" STORE ASSEMBLED MATRICES\n");
				if (nassemblers() == 1) {
					std::string prefix = utils::debugDirectory() + "/assembler";
					std::string suffix = "_aft" + std::to_string(step::ftt.step);
					assembler()->print(builder, prefix.c_str(), suffix.c_str());
				} else {
					for (int i = 0; i < nassemblers(); i++) {
						std::string prefix = utils::debugDirectory() + "/assembler" + std::to_string(i + 1);
						std::string suffix = "_aft" + std::to_string(step::ftt.step);
						assembler()->print(builder, prefix.c_str(), suffix.c_str());
					}
				}
			}
		}

		step::step.type = step::TYPE::FREQUENCY;
	}

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE SOLVER MATRICES\n");
		if (nassemblers() == 1) {
			std::string prefix = utils::debugDirectory() + "/solver";
			solver()->printData(builder, prefix.c_str());
		} else {
			for (int i = 0; i < nsolvers(); i++) {
				std::string prefix = utils::debugDirectory() + "/solver" + std::to_string(i + 1);
				solver()->printData(builder, prefix.c_str());
			}
		}
	}
	eslog::checkpointln("PHYSICS SOLVER: MATRICES ASSEMBLED");
}

void LinearSystem::setDirichlet()
{
	solver()->setDirichlet(builder);
	eslog::checkpointln("PHYSICS SOLVER: DIRICHLET SET");
}

bool LinearSystem::solve()
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE MATRICES FOR LINEAR SOLVER\n");
		std::string prefix = utils::debugDirectory() + "/linsolver";
		solver()->printData(builder, prefix.c_str());
	}

	bool ret = true;

	if (builder->matrices & Builder::Request::KCM) {
		ret = solver()->linearSolver->update();
		if (!ret) return false;
	}
	
	ret = autoopt::solver::evaluate([&] () {
		return solver()->linearSolver->solve();
	});
	if (!ret) return false;


	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE LINEAR SOLVER SOLUTION\n");
		std::string prefix = utils::debugDirectory() + "/linsolver";
		solver()->printSolution(builder, prefix.c_str());
	}
	eslog::checkpointln("PHYSICS SOLVER: SOLUTION COMPUTED");

	return ret;
}

void LinearSystem::solutionChanged()
{
	_builderUpdateSolution();
	eslog::checkpointln("PHYSICS SOLVER: SOLUTION UPDATED");
}

void LinearSystem::processSolution()
{
	assembler()->composer->kernel->processSolution();
	info::mesh->output->updateSolution();
	eslog::checkpointln("PHYSICS SOLVER: SOLUTION PROCESSED");
}

