
#include "ensight.h"

#include "parser/casefile.h"
#include "parser/geometry.h"
#include "parser/variables.h"

#include "config/ecf/input/input.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.hpp"
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

	eslog::info(" ==================================== ENSIGHT GOLD LOADER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == NUMBER OF READERS %69d == \n", MPITools::subset->acrosssize);
	switch (configuration.loader) {
	case InputConfiguration::LOADER::MPI:            eslog::info(" == READER TYPE %75s == \n", "MPI I/O"); break;
	case InputConfiguration::LOADER::MPI_COLLECTIVE: eslog::info(" == READER TYPE %75s == \n", "COLLECTIVE MPI I/O"); break;
	case InputConfiguration::LOADER::POSIX:          eslog::info(" == READER TYPE %75s == \n", "POSIX"); break;
	}
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == CASEFILE %78s == \n", configuration.path.c_str());

	data = new EnsightData();
	data->casefile = new EnsightCasefile(configuration.path);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
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
	eslog::info(" == GEOMETRY %75.2f MB == \n", (double)geofile.totalSize / 1024 / 1024);

	data->geometry = new EnsightGeometry(geofile);
	data->geometry->scan();
	profiler::synccheckpoint("scan");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY SCANNED");

	data->geometry->parse(nodes, elements, regions);
	profiler::synccheckpoint("parse");
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY PARSED");

	data->datafiles = new AsyncFilePack();

	eslog::info(" == VARIABLES %77s == \n", "");
	for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
		if (data->casefile->variables[v].timeset == -1) {
			data->datafiles->add(data->casefile->variables[v].path);
		} else {
			int ts = data->casefile->variables[v].timeset;
			size_t dot = data->casefile->variables[v].path.find_last_of('.') + 1;
			size_t len = data->casefile->variables[v].path.size();
			data->casefile->variables[v].path.replace(dot, len - dot, len - dot, '0');
			for (int t = data->casefile->timesets[ts].tstart; t < data->casefile->timesets[ts].tend; t += data->casefile->timesets[ts].tinc) {
				std::string filepath = data->casefile->variables[v].path;
				std::string suffix = std::to_string(t);
				filepath.replace(filepath.size() - suffix.size(), suffix.size(), suffix);
				data->datafiles->add(filepath);
			}
		}
	}

	data->datafiles->iread();
	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("read_variables");
	eslog::endln("ENSIGHT PARSER: READING VARIABLES STARTED");
}

void InputEnsight::build(Mesh &mesh)
{
	builder::buildOrderedFEM(nodes, elements, regions, mesh);
}

void InputEnsight::variables(Mesh &mesh)
{
	eslog::startln("ENSIGHT VARIABLES: STARTED", "ENSIGHT VARIABLES");
	eslog::info(" ================================== ENSIGHT VARIABLE LOADER =================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	size_t totalSize = 0;
	for (size_t v = 0, fileindex = 0; v < data->casefile->variables.size(); ++v) {
		if (data->casefile->variables[v].timeset == -1) {
			eslog::info(" == %25s: %14s x %d %20s %17.2f MB == \n",
					data->casefile->variables[v].name.c_str(),
					data->casefile->variables[v].type == EnsightCasefile::Variable::Type::NODE ? "NODE" : "ELEMENT",
					data->casefile->variables[v].dimension,
					"",
					(double)data->datafiles->files[fileindex++]->totalSize / 1024 / 1024);
			totalSize += data->datafiles->files[fileindex]->totalSize;
		} else {
			int ts = data->casefile->variables[v].timeset;
			char timesteps[20];
			snprintf(timesteps, 20, "%d->%d:%d", data->casefile->timesets[ts].tstart, data->casefile->timesets[ts].tend, data->casefile->timesets[ts].tinc);
			eslog::info(" == %25s: %14s x %d %20s %17.2f MB == \n",
					data->casefile->variables[v].name.c_str(),
					data->casefile->variables[v].type == EnsightCasefile::Variable::Type::NODE ? "NODE" : "ELEMENT",
					data->casefile->variables[v].dimension,
					timesteps,
					(double)data->datafiles->files[fileindex]->totalSize / 1024 / 1024);
			totalSize += data->datafiles->files[fileindex]->totalSize * (data->casefile->timesets[ts].tend - data->casefile->timesets[ts].tstart) / data->casefile->timesets[ts].tinc;
			fileindex += (data->casefile->timesets[ts].tend - data->casefile->timesets[ts].tstart) / data->casefile->timesets[ts].tinc;
		}
	}

	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == TOTAL VARIABLES SIZE %63.2f MB == \n", (double)totalSize / 1024 / 1024);

	data->datafiles->wait();
	profiler::synccheckpoint("variables_loaded");
	eslog::checkpointln("ENSIGHT VARIABLES: VARIABLES LOADED");

	EnsightVariables variables(*data->casefile, *data->geometry, *data->datafiles);
	variables.scan();
	profiler::synccheckpoint("scan_variables");
	eslog::checkpointln("ENSIGHT VARIABLES: VARIABLES SCANNED");

	variables.parse(mesh);
	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("parse_variables");
	profiler::syncend("ensight");
	eslog::endln("ENSIGHT VARIABLES: VARIABLES PARSED");
}

