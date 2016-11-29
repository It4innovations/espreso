
#include "factory.h"
#include "composition/grid.h"

#include "../../config/description.h"
#include "esmesh.h"

using namespace espreso::input;

void Generator::generate(const ESPRESOGenerator &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.shape) {
	case GENERATOR_SHAPE::GRID:
		Grid::load(mesh, index, size);
		break;
	case GENERATOR_SHAPE::SPHERE:
		ESINFO(GLOBAL_ERROR) << "Implement sphere generator";
		break;
	}
}



