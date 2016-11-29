
#include "factory.h"
#include "esbasis.h"
#include "esinput.h"
#include "../config/description.h"

#include "../../input/loader.h"
#include "../../assembler/assembler.h"

namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
{
	input::Loader::load(configuration, mesh, configuration.env.MPIrank, configuration.env.MPIsize);
	Assembler::compose(configuration, instance, mesh);
}

void Factory::solve(const std::string &outputFile)
{
	instance->init();

	for (size_t i = 0; i < config::solver::TIME_STEPS; i++) {
		instance->solve(_solution);
	}

	instance->finalize();
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


