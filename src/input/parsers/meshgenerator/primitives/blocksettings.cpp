
#include "blocksettings.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "config/ecf/input/block.h"
#include "basis/evaluator/evaluator.h"

using namespace espreso;

BlockSettings::BlockSettings(const BlockGeneratorConfiguration &configuration)
: projection(Expression(configuration.projection_x, { "x", "y", "z" }), Expression(configuration.projection_y, { "x", "y", "z" }), Expression(configuration.projection_z, { "x", "y", "z" })),
  rotation(Expression(configuration.rotation_x, { "x", "y", "z" }), Expression(configuration.rotation_y, { "x", "y", "z" }), Expression(configuration.rotation_z, { "x", "y", "z" }))
{
	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	start = Triple<esint>(configuration.start_x.evaluator->eval() / MeshGenerator::precision, configuration.start_y.evaluator->eval() / MeshGenerator::precision, configuration.start_z.evaluator->eval() / MeshGenerator::precision);
	end   = Triple<esint>(
			(esint)std::round((configuration.start_x.evaluator->eval() + configuration.length_x.evaluator->eval()) / (.1 * MeshGenerator::precision)) / 10,
			(esint)std::round((configuration.start_y.evaluator->eval() + configuration.length_y.evaluator->eval()) / (.1 * MeshGenerator::precision)) / 10,
			(esint)std::round((configuration.start_z.evaluator->eval() + configuration.length_z.evaluator->eval()) / (.1 * MeshGenerator::precision)) / 10);
}

size_t BlockSettings::preferedDomains(const BlockGeneratorConfiguration &configuration)
{
	return configuration.domains_x * configuration.domains_y * configuration.domains_z;
}

