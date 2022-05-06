
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

AcousticComplexLinear::AcousticComplexLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, scheme{configuration.harmonic_solver, 1}, system{}
{

}

void AcousticComplexLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                               HARMONIC COMPLEX == \n");
	eslog::info(" == PHYSICS                                                                        ACOUSTIC == \n");
	eslog::info(" ============================================================================================= \n");

	Variable::list.global.insert(std::make_pair("FREQUENCY", new FrequencyVariable(frequency)));
	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::FREQUENCY);
}

void AcousticComplexLinear::run(step::Step &step)
{
	initSystem(system, this);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
	scheme.init(system);
	assembler.connect(scheme);
	scheme.initFrequency(frequency);

	eslog::info("\n ============================================================================================= \n");
	eslog::info(" = RUN THE SOLVER                       FREQUENCY: MIN %10.4f, MAX %10.4f, STEPS %3d = \n", configuration.harmonic_solver.min_frequency, configuration.harmonic_solver.max_frequency, configuration.harmonic_solver.num_samples);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	system->set(step);
	eslog::info(" ============================================================================================= \n\n");

	while (frequency.current != frequency.final) {
		double start = eslog::time();
		scheme.nextFrequency(frequency);
		eslog::info(" ============================================================================================= \n");
		eslog::info(" = LOAD STEP %2d                                                         FREQUENCY %10.4f = \n", step::step.loadstep + 1, frequency.current);
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

		assembler.evaluate(scheme);
		scheme.composeSystem(frequency, system);
		eslog::info("       = ----------------------------------------------------------------------------- = \n");
		eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

		system->update(step);
		system->solve(step);

		double solution = eslog::time();
		scheme.extractSolution(frequency, system);
		assembler.updateSolution(scheme);
		info::mesh->output->updateSolution(step, frequency);
		eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
		eslog::info("       = ----------------------------------------------------------------------------- = \n");

		eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
	}
}
