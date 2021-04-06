
#include "blocksettings.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "config/ecf/input/block.h"

using namespace espreso;

BlockSettings::BlockSettings(const BlockGeneratorConfiguration &configuration)
: projection(Expression(configuration.projection_x, { "x", "y", "z" }), Expression(configuration.projection_y, { "x", "y", "z" }), Expression(configuration.projection_z, { "x", "y", "z" })),
  rotation(Expression(configuration.rotation_x, { "x", "y", "z" }), Expression(configuration.rotation_y, { "x", "y", "z" }), Expression(configuration.rotation_z, { "x", "y", "z" }))
{
	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	start = Triple<esint>(configuration.start_x / MeshGenerator::precision, configuration.start_y / MeshGenerator::precision, configuration.start_z / MeshGenerator::precision);
	end   = Triple<esint>(
			(configuration.start_x + configuration.length_x) / MeshGenerator::precision,
			(configuration.start_y + configuration.length_y) / MeshGenerator::precision,
			(configuration.start_z + configuration.length_z) / MeshGenerator::precision);
}

size_t BlockSettings::preferedDomains(const BlockGeneratorConfiguration &configuration)
{
	return configuration.domains_x * configuration.domains_y * configuration.domains_z;
}

