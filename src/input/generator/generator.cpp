
#include "composition/grid.h"
#include "composition/gridtower.h"
#include "composition/sphere.h"

#include "generator.h"

#include "../../configuration/input/inputgenerator.h"

using namespace espreso::input;

double Generator::precision = 1e-4;

void Generator::generate(const ESPRESOGenerator &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.shape) {
	case GENERATOR_SHAPE::GRID:
		Grid::load(configuration.grid, mesh, index, size);
		break;
	case GENERATOR_SHAPE::GRID_TOWER:
		GridTower::load(configuration.grid_tower, mesh, index, size);
		break;
	case GENERATOR_SHAPE::SPHERE:
		Sphere::load(configuration.sphere, mesh, index, size);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented generator";
	}
}



