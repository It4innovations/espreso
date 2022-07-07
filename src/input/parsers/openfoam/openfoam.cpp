
#include "openfoam.h"

#include "parser/faceList.h"
#include "parser/labelList.h"
#include "parser/vectorField.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/sysutils.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/builders/builder.h"
#include "wrappers/mpi/communication.h"

#include <fstream>

using namespace espreso;

InputOpenFoam::InputOpenFoam()
: loader(nullptr)
{
	FoamFile::init();
}

InputOpenFoam::~InputOpenFoam()
{
	if (loader) delete loader;
	FoamFile::finish();
}

static int numberOfSubdomains(const std::string &path)
{
	int domains = 1;
	if (info::mpi::rank == 0 && utils::exists(path)) {
		std::ifstream is(path);

		char line[256];
		while (is.getline(line, 256)) {
			if (memcmp(line, "numberOfSubdomains", 18) == 0) {
				domains = std::strtod(line + 18, nullptr);
			}
		}
	}
	Communication::allReduce(&domains, nullptr, 1, MPI_INT, MPI_MAX);
	return domains;
}

static esint numberOfCells(const FoamFileHeader &header)
{
	std::vector<std::string> info = Parser::split(header.note, " ");
	for (size_t i = 0; i < info.size(); ++i) {
		if (StringCompare::caseSensitivePreffix("nCells", info[i])) {
			return std::strtol(Parser::split(info[i], ":")[1].c_str(), nullptr, 10);
		}
	}
	return 0;
}

void InputOpenFoam::load(const InputConfiguration &configuration)
{
	eslog::startln("OPENFOAM PARSER: STARTED", "OPENFOAM PARSER");

	int domains = numberOfSubdomains(configuration.path + "/system/decomposeParDict");
	if (domains == 1) {
		loader = new InputOpenFoamSequential();
	} else {
		loader = new InputOpenFoamParallel(domains);
	}
	loader->load(configuration);

	eslog::endln("OPENFOAM PARSER: PARSED");
}

void InputOpenFoamSequential::load(const InputConfiguration &configuration)
{
	eslog::info(" ================================ SEQUENTIAL OPENFOAM LOADER ================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == MESH %82s == \n", (configuration.path + "/constant/polyMesh/").c_str());

	InputFilePack inputPack;

	OpenFOAMVectorField points(inputPack.add(configuration.path + "/constant/polyMesh/points"));
	OpenFOAMFaceList faces(inputPack.add(configuration.path + "/constant/polyMesh/faces"));

	OpenFOAMLabelList owner(inputPack.add(configuration.path + "/constant/polyMesh/owner"));
	OpenFOAMLabelList neighbour(inputPack.add(configuration.path + "/constant/polyMesh/neighbour"));

	inputPack.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READER PREPARED");

	inputPack.read();
	profiler::synccheckpoint("geometry_read");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READ");

	points.scan(); faces.scan(); owner.scan(); neighbour.scan();
	FoamFile::synchronize({ &points, &faces, &owner, &neighbour });
	profiler::synccheckpoint("scan");

	points.parse(this->nodes.coordinates);
	faces.parse(this->faces.ftype, this->faces.fnodes);
	owner.parse(this->faces.owner);
	neighbour.parse(this->faces.neighbor);

	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");
}

void InputOpenFoamParallel::load(const InputConfiguration &configuration)
{
	eslog::info(" ================================= PARALLEL OPENFOAM LOADER =================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == NUMBER OF SUBDOMAINS %66d == \n", domains);
	eslog::info(" == MESH %82s == \n", (configuration.path + "/processor*/constant/polyMesh/").c_str());

	esint localCells = 0;
	std::vector<esint> nsize(domains / info::mpi::size + 1), esize(domains / info::mpi::size + 1);
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		std::string dict = configuration.path + "/processor" + std::to_string(d) + "/constant/polyMesh/";
		OpenFOAMFaceList::load(dict + "faces", faces.ftype, faces.fnodes, nodes.coordinates.size());
		nsize[i] = OpenFOAMVectorField::load(dict + "points", nodes.coordinates);

		esize[i] = numberOfCells(OpenFOAMLabelList::load(dict + "owner", faces.owner, localCells));
		OpenFOAMLabelList::load(dict + "neighbour", faces.neighbor, localCells);
		faces.neighbor.resize(faces.owner.size(), -1);
		localCells += esize[i];
	}

	std::vector<esint> offset = { (esint)nodes.coordinates.size(), localCells }, sum(2);
	Communication::exscan(sum, offset);
	esint nsum = 0, esum = 0;
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		nodes.blocks.push_back(DatabaseOffset{offset[0] + nsum, nsum, nsize[i]});
		nsum += nsize[i];
		faces.elements.blocks.push_back(DatabaseOffset{offset[1] + esum, esum, esize[i]});
		esum += nsize[i];
	}
	// from local to global adressing
	for (size_t n = 0; n < faces.fnodes.size(); ++n) {
		faces.fnodes[n] += offset[0];
	}
	for (size_t n = 0; n < faces.owner.size(); ++n) {
		faces.owner[n] += offset[1];
	}
	for (size_t n = 0; n < faces.neighbor.size(); ++n) {
		if (faces.neighbor[n] != -1) {
			faces.neighbor[n] += offset[1];
		}
	}

	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");
}

void InputOpenFoam::build(Mesh &mesh)
{
	loader->build(mesh);
}

void InputOpenFoamSequential::build(Mesh &mesh)
{
//	builder::buildOrderedFVM(this->mesh, mesh);
}

void InputOpenFoamParallel::build(Mesh &mesh)
{
	builder::buildDecomposedFVM(nodes, faces, regions, mesh);
}

void InputOpenFoam::variables(Mesh &mesh)
{
	loader->variables(mesh);
}

void InputOpenFoamSequential::variables(Mesh &mesh)
{

}

void InputOpenFoamParallel::variables(Mesh &mesh)
{

}
