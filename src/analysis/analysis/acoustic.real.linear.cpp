
#include "analysis.h"
#include "acoustic.real.linear.h"

#include "analysis/linearsystem/directsystem.h"
#include "analysis/linearsystem/fetisystem.h"
#include "analysis/linearsystem/mklpdsssystem.h"
#include "analysis/linearsystem/multigridsystem.h"

#include "config/ecf/physics/acoustic.h"

#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "output/output.h"

using namespace espreso;

AX_AcousticRealLinear::AX_AcousticRealLinear(AcousticGlobalSettings &gsettings, AcousticLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, scheme{configuration.harmonic_solver, 1}, system{}
{

}

void AX_AcousticRealLinear::init()
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<double>(configuration.feti); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = new AX_MultigridSystem<double>(configuration.hypre); break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<double>(configuration.mklpdss); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = new AX_DirectSystem<double>(configuration.pardiso); break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = new AX_DirectSystem<double>(configuration.superlu); break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = new AX_DirectSystem<double>(configuration.wsmp); break;
	}
	system->init(this);
	scheme.init(system);
	assembler.init(scheme);
}

void AX_AcousticRealLinear::run()
{
	step::Frequency frequency;
	scheme.initFrequency(frequency);

	while (frequency.current != frequency.final) {
		scheme.nextFrequency(frequency);

		assembler.next();
		scheme.composeSystem(frequency, system);

		assembler.fillDirichlet(*system->assembler.dirichlet);
		scheme.composeDirichlet(system);

		scheme.storeScheme(frequency);

		system->update(assembler);
		system->solve();

		scheme.extractSolution(system);
		scheme.storeSolution(frequency);

		assembler.updateSolution();
		info::mesh->output->updateSolution();
	}
}
