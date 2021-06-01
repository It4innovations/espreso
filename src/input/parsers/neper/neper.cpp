
#include "neper.h"
#include "parser/msh.parser.h"
#include "config/ecf/input/input.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

using namespace espreso;

NeperLoader::NeperLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void NeperLoader::load()
{
	eslog::startln("NEPER PARSER: STARTED", "NEPER PARSER");
	profiler::syncstart("neper");

	InputFilePack meshfile;
	meshfile.commitFiles({ _configuration.path });
	meshfile.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("NEPER PARSER: MESH READER PREPARED");

	meshfile.read();
	profiler::synccheckpoint("read");
	eslog::checkpointln("NEPER PARSER: MESH READ");

	meshfile.next();
	DistributedScanner::align(meshfile, "\n");

	NeperMshMesh mesh(meshfile);

	mesh.parse(*this);
	body.resize(etype.size());
	material.resize(etype.size());
	profiler::synccheckpoint("parse");
	profiler::syncend("neper");
	eslog::endln("NEPER PARSER: GEOMETRY PARSED");
}
