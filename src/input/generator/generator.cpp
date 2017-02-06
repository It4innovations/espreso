
#include "composition/grid.h"
#include "composition/sphere.h"

#include "generator.h"
#include "../../configuration/inputgenerator.h"

using namespace espreso::input;

void Generator::generate(const ESPRESOGenerator &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.shape) {
	case GENERATOR_SHAPE::GRID:
		Grid::load(configuration.grid, mesh, index, size);
		break;
	case GENERATOR_SHAPE::SPHERE:
		Sphere::load(configuration.sphere, mesh, index, size);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented generator";
	}
}



