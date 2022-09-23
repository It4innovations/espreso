
#include "meshgenerator.h"
#include "composition/simplegrid.h"
#include "composition/gridgenerator.h"
#include "composition/gridsetgenerator.h"
#include "composition/gridtowergenerator.h"
#include "composition/spheregenerator.h"

#include "basis/logging/profiler.h"
#include "config/ecf/input/generator.h"
#include "esinfo/eslog.h"
#include "input/builders/builder.h"

#include <numeric>

using namespace espreso;

double MeshGenerator::precision = 1e-4;

MeshGenerator::MeshGenerator(InputGeneratorConfiguration &configuration)
: configuration(configuration)
{
	switch (configuration.shape) {
	case INPUT_GENERATOR_SHAPE::GRID:
		generator = new SimpleGridGenerator(configuration.grid);
		break;
	case INPUT_GENERATOR_SHAPE::GRID_SET:
	case INPUT_GENERATOR_SHAPE::GRID_TOWER:
	case INPUT_GENERATOR_SHAPE::SPHERE:
	default:
		eslog::globalerror("Not implemented mesh generator shape.\n");
	}
}

MeshGenerator::~MeshGenerator()
{

}

void MeshGenerator::load(const InputConfiguration &configuration)
{
	profiler::syncstart("mesh_generator");
	generator->nodes(nodes);
	generator->elements(elements);
	generator->neighbors(nodes);
	generator->regions();
	profiler::syncend("mesh_generator");
}

void MeshGenerator::build(Mesh &mesh)
{
	builder::buildDecomposedFEM(nodes, elements, regions, mesh);
}

double MeshGenerator::nextVariables(Mesh &mesh)
{
	return 0;
}



