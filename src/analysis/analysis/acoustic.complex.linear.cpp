
#include "analysis.h"
#include "acoustic.complex.linear.h"

#include "analysis/linearsystem/linearsystem.hpp"
#include "basis/expression/variable.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/acoustic.h"
#include "output/output.h"

using namespace espreso;

AX_AcousticComplexLinear::AX_AcousticComplexLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, scheme{configuration.harmonic_solver, 1}, system{}
{

}

void AX_AcousticComplexLinear::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                               HARMONIC COMPLEX == \n");
	eslog::info(" == PHYSICS                                                                        ACOUSTIC == \n");
	eslog::info(" ============================================================================================= \n");

	initSystem(system, this);
	scheme.init(system);
	assembler.init(scheme);

	Variable::list.global.insert(std::make_pair("FREQUENCY", nullptr));
	info::mesh->output->updateMonitors(step::TYPE::FREQUENCY);
}

void AX_AcousticComplexLinear::run(step::Step &step)
{
	step::Frequency frequency;
	scheme.initFrequency(frequency);
	Variable::list.global["FREQUENCY"] = new FrequencyVariable(frequency);

	eslog::info("\n ============================================================================================= \n");
	eslog::info(" = RUN THE SOLVER                       FREQUENCY: MIN %10.4f, MAX %10.4f, STEPS %3d = \n", configuration.harmonic_solver.min_frequency, configuration.harmonic_solver.max_frequency, configuration.harmonic_solver.num_samples);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	system->info();
	system->set(step);
	eslog::info(" ============================================================================================= \n\n");

	while (frequency.current != frequency.final) {
		double start = eslog::time();
		scheme.nextFrequency(frequency);
		eslog::info(" ============================================================================================= \n");
		eslog::info(" = LOAD STEP %2d                                                         FREQUENCY %10.4f = \n", step::step.loadstep + 1, frequency.current);
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

		assembler.evaluate();
		scheme.composeSystem(frequency, system);
		eslog::info("       = ----------------------------------------------------------------------------- = \n");
		eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

		system->update(step);
		system->solve(step);

		double solution = eslog::time();
		scheme.extractSolution(frequency, system);
		assembler.updateSolution();
		info::mesh->output->updateSolution(step, frequency);
		eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
		eslog::info("       = ----------------------------------------------------------------------------- = \n");

		eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
	}
}
