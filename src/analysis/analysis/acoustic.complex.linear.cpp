
#include "analysis.h"
#include "acoustic.complex.linear.h"

#include "analysis/linearsystem/directsystem.h"
#include "analysis/linearsystem/fetisystem.h"
#include "analysis/linearsystem/mklpdsssystem.h"
#include "analysis/linearsystem/multigridsystem.h"

#include "basis/expression/variable.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/acoustic.h"
#include "output/output.h"

using namespace espreso;

AX_AcousticComplexLinear::AX_AcousticComplexLinear(AcousticGlobalSettings &gsettings, AcousticLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, scheme{configuration.harmonic_solver, 1}, system{}
{

}

void AX_AcousticComplexLinear::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                               HARMONIC COMPLEX == \n");
	eslog::info(" == PHYSICS                                                                        ACOUSTIC == \n");
	eslog::info(" ============================================================================================= \n");
//	switch (configuration.solver) {
//	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<std::complex<double>>(configuration.feti); break;
//	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = new AX_MultigridSystem<std::complex<double>>(configuration.hypre); break;
//	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<std::complex<double>>(configuration.mklpdss); break;
//	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = new AX_DirectSystem<std::complex<double>>(configuration.pardiso); break;
//	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = new AX_DirectSystem<std::complex<double>>(configuration.superlu); break;
//	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = new AX_DirectSystem<std::complex<double>>(configuration.wsmp); break;
//	}
	system = new AX_MKLPDSSSystem<AX_AcousticComplexLinear>(this, configuration.mklpdss);
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
