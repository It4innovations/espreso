
#include "analysis.h"
#include "acoustic.complex.linear.h"

#include "analysis/linearsystem/directsystem.h"
#include "analysis/linearsystem/fetisystem.h"
#include "analysis/linearsystem/mklpdsssystem.h"
#include "analysis/linearsystem/multigridsystem.h"

#include "basis/expression/variable.h"
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

	while (frequency.current != frequency.final) {
		scheme.nextFrequency(frequency);

		assembler.evaluate();
		scheme.composeSystem(frequency, system);

		system->info();
		system->set(step);
		system->update(step);
		system->solve(step);

		scheme.extractSolution(frequency, system);

		assembler.updateSolution();
		info::mesh->output->updateSolution(step, frequency);
	}
}
