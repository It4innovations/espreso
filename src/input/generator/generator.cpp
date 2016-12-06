
#include "composition/grid.h"
#include "composition/sphere.h"

#include "../../config/description.h"
#include "esmesh.h"
#include "generator.h"

using namespace espreso::input;

void Generator::generate(const GlobalConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.generator.shape) {
	case GENERATOR_SHAPE::GRID:
		Grid::load(configuration, mesh, index, size);
		break;
	case GENERATOR_SHAPE::SPHERE:
		Sphere::load(configuration, mesh, index, size);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented generator";
	}
}



