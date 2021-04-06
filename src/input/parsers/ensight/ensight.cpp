
#include "ensight.h"
#include "parser/casefile.h"
#include "parser/geometry.h"

#include "config/ecf/input/input.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "input/meshbuilder.h"

using namespace espreso;

EnsightLoader::EnsightLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void EnsightLoader::load()
{
	eslog::startln("ENSIGHT PARSER: STARTED", "ENSIGHT PARSER");
	profiler::syncstart("ensight");

	EnsightCasefile casefile(_configuration.path);
	profiler::synccheckpoint("casefile");
	eslog::checkpointln("ENSIGHT PARSER: CASEFILE READ");

	InputFilePack geofile;
	geofile.commitFiles({ casefile.geometry });
	geofile.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READER PREPARED");

	geofile.read();
	profiler::synccheckpoint("read");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READ");

	geofile.next();

	EnsightGeometry geometry(geofile);
	geometry.scan();
	profiler::synccheckpoint("scan");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY SCANNED");

	geometry.parse(*this);
	removeDuplicates = true;
	body.resize(etype.size());
	material.resize(etype.size());
	profiler::synccheckpoint("parse");
	profiler::syncend("ensight");
	eslog::endln("ENSIGHT PARSER: GEOMETRY PARSED");
}
