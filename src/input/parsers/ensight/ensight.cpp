
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

namespace espreso {

struct EnsightData {
	EnsightCasefile *casefile;
	EnsightGeometry *geometry;
	AsyncFilePack *datafiles;

	~EnsightData()
	{
		delete casefile;
		delete geometry;
		delete datafiles;
	}
};

}

using namespace espreso;

InputEnsight::~InputEnsight()
{
	delete data;
}

void InputEnsight::load(const InputConfiguration &configuration)
{
	eslog::startln("ENSIGHT PARSER: STARTED", "ENSIGHT PARSER");
	profiler::syncstart("ensight");

	data = new EnsightData();
	data->casefile = new EnsightCasefile(configuration.path);
	profiler::synccheckpoint("casefile");
	eslog::checkpointln("ENSIGHT PARSER: CASEFILE READ");

	InputFilePack geofile({ data->casefile->geometry });
	geofile.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READER PREPARED");

	geofile.read();
	profiler::synccheckpoint("read");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READ");

	geofile.next();

	data->geometry = new EnsightGeometry(geofile);
	data->geometry->scan();
	profiler::synccheckpoint("scan");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY SCANNED");

	data->geometry->parse(mesh);
	profiler::synccheckpoint("parse");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY PARSED");

	data->datafiles = new AsyncFilePack();

	for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
		if (data->casefile->variables[v].time == -1) {
			data->datafiles->add(data->casefile->variables[v].path);
		} else {
			size_t dot = data->casefile->variables[v].path.find_last_of('.') + 1;
			size_t len = data->casefile->variables[v].path.size();
			data->casefile->variables[v].path.replace(dot, len - dot, len - dot, '0');
			for (size_t t = 0; t < data->casefile->times[data->casefile->variables[v].time].size(); ++t) {
				std::string filepath = data->casefile->variables[v].path;
				std::string suffix = std::to_string(t + 1);
				filepath.replace(filepath.size() - suffix.size(), suffix.size(), suffix);
				data->datafiles->add(filepath);
			}
		}
	}

	data->datafiles->iread();
	profiler::synccheckpoint("read_variables");
	eslog::endln("ENSIGHT PARSER: READING VARIABLES STARTED");
}

void InputEnsight::build(Mesh &mesh)
{
	builder::buildOrderedFEM(this->mesh, mesh);
}

void InputEnsight::variables(Mesh &mesh)
{
	eslog::startln("ENSIGHT VARIABLES: STARTED", "ENSIGHT VARIABLES");

	data->datafiles->wait();
	profiler::synccheckpoint("variables_loaded");
	eslog::checkpointln("ENSIGHT VARIABLES: VARIABLES LOADED");

	EnsightVariables variables(*data->casefile, *data->geometry, *data->datafiles);
	variables.scan();
	profiler::synccheckpoint("scan_variables");
	eslog::checkpointln("ENSIGHT VARIABLES: VARIABLES SCANNED");

	variables.parse(mesh);
	profiler::synccheckpoint("parse_variables");
	profiler::syncend("ensight");
	eslog::endln("ENSIGHT VARIABLES: VARIABLES PARSED");
}

