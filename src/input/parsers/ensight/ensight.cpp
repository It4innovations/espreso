
#include "ensight.h"
#include "parser/casefile.h"
#include "parser/geometry.h"
#include "parser/variables.h"

#include "config/ecf/input/input.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "input/builders/builder.h"

using namespace espreso;

void InputEnsight::load(const InputConfiguration &configuration)
{
	eslog::startln("ENSIGHT PARSER: STARTED", "ENSIGHT PARSER");
	profiler::syncstart("ensight");

	EnsightCasefile casefile(configuration.path);
	profiler::synccheckpoint("casefile");
	eslog::checkpointln("ENSIGHT PARSER: CASEFILE READ");

	InputFilePack geofile({ casefile.geometry });
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

	geometry.parse(mesh);
	profiler::synccheckpoint("parse");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY PARSED");

	std::vector<std::string> files;
	for (size_t v = 0; v < casefile.variables.size(); ++v) {
		if (casefile.variables[v].time == -1) {
			files.push_back(casefile.variables[v].path);
		} else {
			size_t dot = casefile.variables[v].path.find_last_of('.') + 1;
			size_t len = casefile.variables[v].path.size();
			casefile.variables[v].path.replace(dot, len - dot, len - dot, '0');
			for (size_t t = 0; t < casefile.times[casefile.variables[v].time].size(); ++t) {
				std::string file = casefile.variables[v].path;
				std::string suffix = std::to_string(t + 1);
				file.replace(file.size() - suffix.size(), suffix.size(), suffix);
				files.push_back(file);
			}
		}
	}

	InputFilePack datafiles(files);

	datafiles.prepare();
	profiler::synccheckpoint("prepare_variables");
	eslog::checkpointln("ENSIGHT PARSER: VARIABLES READER PREPARED");

	datafiles.read();
	profiler::synccheckpoint("read_variables");
	eslog::checkpointln("ENSIGHT PARSER: VARIABLES READ");

	EnsightVariables variables(casefile, geometry, datafiles);
	variables.scan();
	profiler::synccheckpoint("scan_variables");
	eslog::checkpointln("ENSIGHT PARSER: VARIABLES SCANNED");

	variables.parse();
	profiler::synccheckpoint("parse_variables");
	profiler::syncend("ensight");
	eslog::endln("ENSIGHT PARSER: VARIABLES PARSED");
}

void InputEnsight::build(Mesh &mesh)
{
	builder::buildOrderedFEM(this->mesh, mesh);
}

void InputEnsight::variables(Mesh &mesh)
{

}

