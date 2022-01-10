
#include "netgen.h"
#include "parser/neutralmesh.h"
#include "config/ecf/input/input.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

using namespace espreso;

NetgenNeutralLoader::NetgenNeutralLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{


}

void NetgenNeutralLoader::load()
{
	eslog::startln("NETGEN PARSER: STARTED", "NETGEN PARSER");
	profiler::syncstart("netgen");

	InputFilePack meshfile({ _configuration.path });
	meshfile.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("NETGEN PARSER: MESH READER PREPARED");

	meshfile.read();
	profiler::synccheckpoint("read");
	eslog::checkpointln("NETGEN PARSER: MESH READ");

	meshfile.next();
	DistributedScanner::align(meshfile, "\n");

	NetgenNeutralMesh mesh(meshfile);

	mesh.parse(*this);
	body.resize(etype.size());
	material.resize(etype.size());
	profiler::synccheckpoint("parse");
	profiler::syncend("netgen");
	eslog::endln("NETGEN PARSER: GEOMETRY PARSED");
}
