
#include "meshgenerator.h"

#include "composition/gridgenerator.h"
#include "composition/gridtowergenerator.h"
#include "composition/spheregenerator.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/generator.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "mesh/mesh.h"

using namespace espreso;

double MeshGenerator::precision = 1e-4;

MeshGenerator::MeshGenerator(const InputGeneratorConfiguration &configuration)
: MeshBuilder(MeshBuilder::TYPE::GENERATED),
  _configuration(configuration)
{

}

void MeshGenerator::load()
{
	profiler::syncstart("mesh_generator");
	switch (_configuration.shape) {
	case INPUT_GENERATOR_SHAPE::GRID:
		GridGenerator::generate(_configuration.grid, *this);
		info::mesh->preferedDomains = BlockSettings::preferedDomains(_configuration.grid);
		break;
	case INPUT_GENERATOR_SHAPE::GRID_TOWER:
		GridTowerGenerator::generate(_configuration.grid_tower, *this);
		info::mesh->preferedDomains = BlockSettings::preferedDomains(_configuration.grid_tower.grids.at(GridTowerGenerator::gridIndex(_configuration.grid_tower)));
		break;
	case INPUT_GENERATOR_SHAPE::SPHERE:
		SphereGenerator::generate(_configuration.sphere, *this);
		info::mesh->preferedDomains = BlockSettings::preferedDomains(_configuration.sphere);
		break;
	default:
		eslog::globalerror("Not implemented mesh generator shape.\n");
	}
	profiler::syncend("mesh_generator");
}



