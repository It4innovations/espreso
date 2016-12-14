
#include "factory.h"
#include "esbasis.h"
#include "esinput.h"
#include "../../input/loader.h"
#include "../../assembler/assembler.h"
#include "../../config/globalconfiguration.h"

namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
{
	input::Loader::load(configuration, mesh, configuration.env.MPIrank, configuration.env.MPIsize);
	Assembler::compose(configuration, instance, mesh);
}

void Factory::solve(const std::string &outputFile)
{
	instance->init();
	instance->solve(_solution);
	instance->finalize();
}

void Factory::check(const Results &configuration)
{
	auto norm = [&] () {
		double n = 0, sum = 0;
		for (size_t i = 0; i < _solution.size(); i++) {
			for (size_t j = 0; j < _solution[i].size(); j++) {
				n += _solution[i][j] * _solution[i][j];
			}
		}

		MPI_Allreduce(&n, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return sqrt(sum);
	};

	if (configuration.norm != 0) {
		double nn = norm();
		ESTEST(EVALUATION)
			<< (fabs(nn - configuration.norm) > 1e-3 && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
			<< "Norm of the solution " << nn << " is not " << configuration.norm << ".";
	}
}

double Factory::norm() const
{
	double n = 0, sum = 0;
	for (size_t i = 0; i < _solution.size(); i++) {
		for (size_t j = 0; j < _solution[i].size(); j++) {
			n += _solution[i][j] * _solution[i][j];
		}
	}

	MPI_Allreduce(&n, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(sum);
}

}


