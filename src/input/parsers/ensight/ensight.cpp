
#include "ensight.h"

#include "parser/casefile.h"
#include "parser/geometry.h"
#include "parser/variables.h"

#include "config/ecf/input/input.h"
#include "basis/io/inputfile.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "esinfo/stepinfo.h"
#include "input/builders/builder.h"

namespace espreso {

struct EnsightData {
	EnsightCasefile *casefile;
	EnsightGeometry *geometry;
	AsyncFilePack *datafiles;
	EnsightVariables *variables;

	~EnsightData()
	{
		delete casefile;
		delete geometry;
		delete datafiles;
		delete variables;
	}
};

}

using namespace espreso;

InputEnsight::~InputEnsight()
{
	delete data;
}

InputEnsight::InputEnsight()
: data(nullptr), timeset(0), timestep(0)
{

}

void InputEnsight::load(const InputConfiguration &configuration)
{
	eslog::startln("ENSIGHT PARSER: STARTED", "ENSIGHT PARSER");

	eslog::info(" ==================================== ENSIGHT GOLD LOADER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == CASEFILE %78s == \n", configuration.path.c_str());

	data = new EnsightData();
	data->casefile = new EnsightCasefile(configuration.path);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::checkpointln("ENSIGHT PARSER: CASEFILE READ");

	InputFilePack geofile({ data->casefile->geometry });
	geofile.prepare();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READER PREPARED");

	geofile.read();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READ");

	geofile.next();
	eslog::info(" == GEOMETRY %75.2f MB == \n", (double)geofile.totalSize / 1024 / 1024);

	data->geometry = new EnsightGeometry(geofile);
	data->geometry->scan();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY SCANNED");

	data->geometry->parse(nodes, elements, regions);
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY PARSED");

	data->datafiles = new AsyncFilePack();

	if (configuration.ensight.variables.size()) {
		auto variableNames = Parser::split(configuration.ensight.variables, ",");
		for (size_t v = 0; v < variableNames.size(); ++v) {
			variableNames[v] = Parser::strip(variableNames[v]);
		}
		for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
			data->casefile->variables[v].skip = true;
		}
		for (size_t i = 0; i < variableNames.size(); ++i) {
			for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
				if (data->casefile->variables[v].skip && StringCompare::caseInsensitiveEq(variableNames[i], data->casefile->variables[v].name)) {
					data->casefile->variables[v].skip = false;
				}
			}
		}
	}

	data->variables = new EnsightVariables(*data->casefile, *data->geometry, *data->datafiles);

	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	timeset = configuration.ensight.timeset - 1;
	eslog::info(" == TIME SET         %70d == \n", timeset + 1);
	eslog::info(" == TIME STEPS       %70d == \n", (data->casefile->timesets[timeset].tend - data->casefile->timesets[timeset].tstart) / data->casefile->timesets[timeset].tinc + 1);
	eslog::info(" == > START          %70g == \n", data->casefile->timesets[timeset].values.front());
	eslog::info(" == > LAST           %70g == \n", data->casefile->timesets[timeset].values.back());
	eslog::info(" == VARIABLES %77s == \n", "");
	for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
		if (!data->casefile->variables[v].skip) {
			switch (data->casefile->variables[v].type) {
			case EnsightCasefile::Variable::Type::ELEMENT:
				eslog::info(" == > %s %*s x %d == \n", data->casefile->variables[v].name.c_str(), 80 - data->casefile->variables[v].name.size(), "NODE", data->casefile->variables[v].dimension);
				variables.elements.push_back({ data->casefile->variables[v].dimension, data->casefile->variables[v].name });
				break;
			case EnsightCasefile::Variable::Type::NODE:
				eslog::info(" == > %s %*s x %d == \n", data->casefile->variables[v].name.c_str(), 80 - data->casefile->variables[v].name.size(), "ELEMENT", data->casefile->variables[v].dimension);
				variables.nodes.push_back({ data->casefile->variables[v].dimension, data->casefile->variables[v].name });
				break;
			}
		}
	}

	for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
		if (!data->casefile->variables[v].skip) {
			if (data->casefile->variables[v].timeset == -1) {
				data->datafiles->add(data->casefile->variables[v].path);
			}
			if (data->casefile->variables[v].timeset == timeset) {
				size_t dot = data->casefile->variables[v].path.find_last_of('.') + 1;
				size_t len = data->casefile->variables[v].path.size();
				data->casefile->variables[v].path.replace(dot, len - dot, len - dot, '0');
				std::string filepath = data->casefile->variables[v].path;
				std::string suffix = std::to_string(data->casefile->timesets[timeset].tstart);
				filepath.replace(filepath.size() - suffix.size(), suffix.size(), suffix);
				data->datafiles->add(filepath);
			}
		}
	}

	data->datafiles->iread(MPITools::subset->within);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("ENSIGHT PARSER: READING VARIABLES STARTED");
}

void InputEnsight::build(Mesh &mesh)
{
	builder::buildOrderedFEM(nodes, elements, regions, mesh);
}

void InputEnsight::initVariables(Mesh &mesh)
{
	eslog::startln("ENSIGHT VARIABLES LOADER: STARTED", "VARIABLES LOADER");
}

void InputEnsight::finishVariables()
{
	builder::orderedValuesFinish(variables);
	eslog::endln("VARIABLES LOADER: FINISHED");
}

int InputEnsight::timeSteps()
{
	return (data->casefile->timesets[timeset].tend - data->casefile->timesets[timeset].tstart) / data->casefile->timesets[timeset].tinc + 1;
}

void InputEnsight::nextTimeStep(Mesh &mesh)
{
	eslog::start("TIME STEP LOADER: STARTED", "TIME STEP LOADER");
	step::time.current = data->casefile->timesets[timeset].values[timestep];
	eslog::info(" == TIME STEP %77g == \n", step::time.current);

	data->datafiles->wait(MPITools::subset->within);
	eslog::checkpointln("TIME STEP LOADER: VARIABLES LOADED");

	if (timestep == 0) {
		data->variables->scan();
		eslog::checkpointln("VARIABLES LOADER: VARIABLES SCANNED");
	}

	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		variables.nodes[v].data.clear();
	}
	for (size_t v = 0; v < variables.elements.size(); ++v) {
		variables.elements[v].data.clear();
	}
	data->variables->parse(mesh, variables);
	eslog::checkpointln("TIME STEP LOADER: VARIABLES PARSED");

	if (timestep == 0) {
		builder::orderedValuesInit(variables, mesh);
	}

	if (++timestep < data->casefile->timesets[timeset].values.size()) {
		data->datafiles->clear();
		for (size_t v = 0; v < data->casefile->variables.size(); ++v) {
			if (!data->casefile->variables[v].skip) {
				if (data->casefile->variables[v].timeset == timeset) {
					size_t dot = data->casefile->variables[v].path.find_last_of('.') + 1;
					size_t len = data->casefile->variables[v].path.size();
					data->casefile->variables[v].path.replace(dot, len - dot, len - dot, '0');
					std::string filepath = data->casefile->variables[v].path;
					std::string suffix = std::to_string(data->casefile->timesets[timeset].tstart + data->casefile->timesets[timeset].tinc * timestep);
					filepath.replace(filepath.size() - suffix.size(), suffix.size(), suffix);
					data->datafiles->add(filepath);
				}
			}
		}
		data->datafiles->iread(MPITools::subset->within);
		eslog::checkpointln("TIME STEP LOADER: IREAD NEXT TIME STEP");
	}

	builder::orderedValuesNext(variables, mesh);
	eslog::endln("TIME STEP LOADER: FINISHED");
}

