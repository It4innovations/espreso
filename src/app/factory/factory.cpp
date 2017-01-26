
#include "factory.h"
#include "esbasis.h"
#include "esinput.h"
#include "../../input/loader.h"
#include "../../assembler/assembler.h"
#include "../../config/globalconfiguration.h"
#include "../../mesh/settings/evaluator.h"

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
	double epsilon = 1e-2;
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
			<< (fabs(nn - configuration.norm) > 1e-2 && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
			<< (fabs(nn - configuration.norm) > epsilon && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
			<< "Norm of the solution " << nn << " is not " << configuration.norm << ".";
	}

	auto evaluateProperty = [&] (const std::string &value, Property property, size_t DOF) {
		if (instance->physics().pointDOFs[DOF] == property && value.size()) {
			CoordinatesEvaluator evaluator(value, mesh.coordinates());
			for (size_t p = 0; p < mesh.parts(); p++) {
				for (size_t n = 0; n < mesh.coordinates().localSize(p); n++) {
					eslocal index = mesh.coordinates().localToCluster(p)[n];
					ESTEST(EVALUATION)
						<< (fabs(evaluator.evaluate(mesh.coordinates()[index]) - _solution[p][mesh.nodes()[index]->DOFIndex(p, DOF)]) > epsilon ? TEST_FAILED : TEST_PASSED)
						<< "Incorrect " << property << " of the solution.";
				}
			}
		}
	};

	for (size_t DOF = 0; DOF < instance->physics().pointDOFs.size(); DOF++) {
		evaluateProperty(configuration.displacement_x, Property::DISPLACEMENT_X, DOF);
		evaluateProperty(configuration.displacement_y, Property::DISPLACEMENT_Y, DOF);
		evaluateProperty(configuration.displacement_z, Property::DISPLACEMENT_Z, DOF);
		evaluateProperty(configuration.temperature   , Property::TEMPERATURE   , DOF);
		evaluateProperty(configuration.pressure      , Property::PRESSURE      , DOF);
	};
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


