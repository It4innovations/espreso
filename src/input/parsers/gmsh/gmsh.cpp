
#include "gmsh.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "wrappers/gmsh/w.gmsh.h"

using namespace espreso;

GMSHGenerator::GMSHGenerator(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void GMSHGenerator::load()
{
	eslog::startln("GMSH GENERATOR: STARTED", "GMSH GENERATOR");
	profiler::syncstart("gmsh");

	GMSH::generate(*this);

	body.resize(etype.size());
	material.resize(etype.size());
	profiler::syncend("gmsh");
	eslog::endln("GMSH GENERATOR: MESH GENERATED");
}

