
#include "blocksettings.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "config/ecf/input/block.h"
#include "basis/evaluator/evaluator.h"

using namespace espreso;

BlockSettings::BlockSettings(const BlockGeneratorConfiguration &configuration)
{
	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	start = Triple<esint>(configuration.start_x.evaluator->evaluate() / MeshGenerator::precision, configuration.start_y.evaluator->evaluate() / MeshGenerator::precision, configuration.start_z.evaluator->evaluate() / MeshGenerator::precision);
	end   = Triple<esint>(
			(esint)std::round((configuration.start_x.evaluator->evaluate() + configuration.length_x.evaluator->evaluate()) / (.1 * MeshGenerator::precision)) / 10,
			(esint)std::round((configuration.start_y.evaluator->evaluate() + configuration.length_y.evaluator->evaluate()) / (.1 * MeshGenerator::precision)) / 10,
			(esint)std::round((configuration.start_z.evaluator->evaluate() + configuration.length_z.evaluator->evaluate()) / (.1 * MeshGenerator::precision)) / 10);

	projection_x = configuration.projection_x.evaluator;
	projection_y = configuration.projection_y.evaluator;
	projection_z = configuration.projection_z.evaluator;
}

size_t BlockSettings::preferedDomains(const BlockGeneratorConfiguration &configuration)
{
	return configuration.domains_x * configuration.domains_y * configuration.domains_z;
}

