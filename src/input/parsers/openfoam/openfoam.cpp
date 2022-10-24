
#include "openfoam.h"

#include "parser/faceList.h"
#include "parser/labelList.h"
#include "parser/vectorField.h"
#include "parser/polyBoundaryMesh.h"
#include "parser/time.h"

#include "basis/containers/tarray.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "input/builders/builder.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

#include <fstream>
#include <numeric>
#include <algorithm>
#include <sstream>

using namespace espreso;

InputOpenFoam::InputOpenFoam()
: timestep(0), loader(nullptr)
{

}

InputOpenFoam::~InputOpenFoam()
{
	if (loader) delete loader;
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
	ParallelFoamFile::init();

	int domains = numberOfSubdomains(configuration.path + "/system/decomposeParDict");
	if (domains == 1) {
		loader = new InputOpenFoamSequential();
	} else {
		if (configuration.openfoam.direct_load) {
			loader = new InputOpenFoamParallelDirect(domains);
		} else {
			loader = new InputOpenFoamParallel(domains);
		}
	}
	loader->load(configuration);

	ParallelFoamFile::finish();
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
	ParallelFoamFile::synchronize({ &points, &faces, &owner, &neighbour });
	profiler::synccheckpoint("scan");

	this->faces.edist.blocks.push_back(DatabaseOffset{0, 0, numberOfCells(owner.header)});

	points.parse(this->nodes.coordinates);
	faces.parse(this->faces.ftype, this->faces.fnodes);
	owner.parse(this->faces.owner);
	neighbour.parse(this->faces.neighbor);

	ivariables(configuration);
	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");
}

void InputOpenFoamParallel::load(const InputConfiguration &configuration)
{
	eslog::info(" ================================ PARALLEL OPENFOAM LOADER ==================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == NUMBER OF SUBDOMAINS %66d == \n", domains);
	eslog::info(" == MESH %82s == \n", (configuration.path + "/processor*/constant/polyMesh/").c_str());

	faces.elements = 0;
	esint nsum = 0, esum = 0;
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		std::string dict = configuration.path + "/processor" + std::to_string(d) + "/constant/polyMesh/";
		OpenFOAMFaceList::load(dict + "faces", faces.ftype, faces.fnodes, nodes.coordinates.size());
		nsum += OpenFOAMVectorField::load(dict + "points", nodes.coordinates);

		esum += numberOfCells(OpenFOAMLabelList::load(dict + "owner", faces.owner, faces.elements));
		OpenFOAMLabelList::load(dict + "neighbour", faces.neighbor, faces.elements);
		faces.neighbor.resize(faces.owner.size(), -1);
		faces.elements = esum;
	}

	esint offset[2] = { nsum, esum };
	Communication::exscan<esint>(offset, nullptr, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	nodes.blocks.push_back(DatabaseOffset{offset[0], 0, nsum});
	faces.edist.blocks.push_back(DatabaseOffset{offset[1], 0, esum});
	variables.ndist.blocks.push_back(nodes.blocks.back());
	variables.edist.blocks.push_back(faces.edist.blocks.back());

	// from local to global addressing
	for (size_t f = 0, foffset = 0; f < faces.ftype.size(); foffset += Element::encode(faces.ftype[f++]).nodes) {
		PolyElement poly(faces.ftype[f], faces.fnodes.data() + foffset);
		for (int n = 0; n < Element::encode(faces.ftype[f]).nodes; ++n) {
			if (poly.isNode(n)) {
				faces.fnodes[n + foffset] += offset[0];
			}
		}
	}
	for (size_t n = 0; n < faces.owner.size(); ++n) {
		faces.owner[n] += offset[1];
	}
	for (size_t n = 0; n < faces.neighbor.size(); ++n) {
		if (faces.neighbor[n] != -1) {
			faces.neighbor[n] += offset[1];
		}
	}

	ivariables(configuration);
	eslog::info(" ============================================================================================= \n\n");
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");
}

void InputOpenFoamParallelDirect::load(const InputConfiguration &configuration)
{
	eslog::info(" == NUMBER OF SUBDOMAINS %66d == \n", domains);
	eslog::info(" == MESH %82s == \n", (configuration.path + "/processor*/constant/polyMesh/").c_str());

	faces.elements = 0;
	std::vector<esint> nsize(domains / info::mpi::size + 1), esize(domains / info::mpi::size + 1);
	std::vector<OpenFOAMPolyBoundaryMesh::OpenFOAMBoundary> neighbors;
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		std::string dict = configuration.path + "/processor" + std::to_string(d) + "/constant/polyMesh/";

		OpenFOAMPolyBoundaryMesh regions(dict + "boundary");
		OpenFOAMFaceList::load(dict + "faces", faces.ftype, faces.fnodes, nodes.coordinates.size());
		OpenFOAMLabelList::load(dict + "pointProcAddressing", nodes.ids, nodes.coordinates.size());
		nsize[i] = OpenFOAMVectorField::load(dict + "points", nodes.coordinates);

		esize[i] = numberOfCells(OpenFOAMLabelList::load(dict + "owner", faces.owner, faces.elements));
		OpenFOAMLabelList::load(dict + "neighbour", faces.neighbor, faces.elements);
		faces.neighbor.resize(faces.owner.size(), -1);

		for (size_t b = 0; b < regions.boundaries.size(); ++b) {
			if (regions.boundaries[b].type == OpenFOAMPolyBoundaryMesh::OpenFOAMBoundaryType::processor) {
				neighbors.push_back(regions.boundaries[b]);
				neighbors.back().startFace += faces.elements;
			}
		}
		faces.elements += esize[i];
	}

	std::sort(neighbors.begin(), neighbors.end(), [] (const OpenFOAMPolyBoundaryMesh::OpenFOAMBoundary &b1, const OpenFOAMPolyBoundaryMesh::OpenFOAMBoundary &b2) { return b1.neighbor < b2.neighbor; });

	ivector<esint> fdist(faces.ftype.size() + 1);
	fdist[0] = 0;
	for (size_t e = 0; e < faces.ftype.size(); ++e) {
		fdist[e + 1] = fdist[e] + Element::encode(faces.ftype[e]).nodes;
	}

	nodes.ranks.distribution.resize(nodes.coordinates.size() + 2, 1);
	nodes.ranks.distribution.front() = 0;
	std::vector<int> lastRank(nodes.coordinates.size(), info::mpi::rank);
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		for (size_t n = 0; n < neighbors.size(); ++n) {
			for (esint f = neighbors[n].startFace; f < neighbors[n].startFace + neighbors[n].nFaces; ++f) {
				esint noffset = Element::encode(faces.ftype[f]).code == Element::CODE::POLYGON ? 1 : 0;
				for (esint node = fdist[f] + noffset; node < fdist[f + 1]; ++node) {
					if (lastRank[faces.fnodes[node]] != neighbors[n].neighbor) {
						lastRank[faces.fnodes[node]] = neighbors[n].neighbor;
						++nodes.ranks.distribution[faces.fnodes[node] + 1];
					}
				}
			}
		}
	}

	utils::sizesToOffsets(nodes.ranks.distribution.data(), nodes.ranks.distribution.data() + nodes.ranks.distribution.size());
	nodes.ranks.data.resize(nodes.ranks.distribution.back(), info::mpi::rank);
	std::fill(lastRank.begin(), lastRank.end(), info::mpi::rank);
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		size_t n = 0;
		auto insert = [&] (const OpenFOAMPolyBoundaryMesh::OpenFOAMBoundary &neighbor) {
			for (esint f = neighbor.startFace; f < neighbor.startFace + neighbor.nFaces; ++f) {
				esint noffset = Element::encode(faces.ftype[f]).code == Element::CODE::POLYGON ? 1 : 0;
				for (esint node = fdist[f] + noffset; node < fdist[f + 1]; ++node) {
					if (lastRank[faces.fnodes[node]] != neighbors[n].neighbor) {
						lastRank[faces.fnodes[node]] = neighbors[n].neighbor;
						nodes.ranks.data[nodes.ranks.distribution[faces.fnodes[node] + 1]++] = neighbor.neighbor;
					}
				}
			}
		};
		for (; n < neighbors.size() && neighbors[n].neighbor < info::mpi::rank; ++n) {
			insert(neighbors[n]);
		}
		for (size_t c = 0; c < nodes.coordinates.size(); ++c) {
			++nodes.ranks.distribution[c + 1];
		}
		for (; n < neighbors.size(); ++n) {
			insert(neighbors[n]);
		}
	}
	nodes.ranks.distribution.pop_back();

	for (size_t n = 0; n < neighbors.size(); ++n) {
		if (nodes.neighbors.size() == 0 || nodes.neighbors.back() != neighbors[n].neighbor) {
			if (neighbors[n].neighbor != info::mpi::rank) {
				nodes.neighbors.push_back(neighbors[n].neighbor);
			}
		}
	}

//	ivariables(configuration);
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");
}

void InputOpenFoam::build(Mesh &mesh)
{
	loader->build(mesh);
}

void InputOpenFoamSequential::build(Mesh &mesh)
{
	builder::buildOrderedFVM(nodes, faces, regions, mesh);
}

void InputOpenFoamParallel::build(Mesh &mesh)
{
	builder::buildChunkedFVM(nodes, faces, regions, mesh);
}

void InputOpenFoamParallelDirect::build(Mesh &mesh)
{
	builder::buildDecomposedFVM(nodes, faces, regions, mesh);
}

void InputOpenFoam::initVariables(Mesh &mesh)
{
	loader->initVariables(mesh);
}

void InputOpenFoamSequential::initVariables(Mesh &mesh)
{

}

void InputOpenFoamParallel::initVariables(Mesh &mesh)
{
	eslog::startln("OPENFOAM VARIABLES LOADER: STARTED", "VARIABLES LOADER");
	builder::chunkedValuesInit(variables, mesh);
}

void InputOpenFoamParallelDirect::initVariables(Mesh &mesh)
{

}

void InputOpenFoam::finishVariables()
{
	loader->finishVariables();
}

void InputOpenFoamSequential::finishVariables()
{

}

void InputOpenFoamParallel::finishVariables()
{
	builder::chunkedValuesFinish(variables);
	eslog::endln("VARIABLES LOADER: FINISHED");
}

void InputOpenFoamParallelDirect::finishVariables()
{

}

int InputOpenFoam::variables()
{
	return loader->timesteps.size();
}

void InputOpenFoam::nextVariables(Mesh &mesh)
{
	return loader->nextVariables(mesh);
}

void InputOpenFoamSequential::nextVariables(Mesh &mesh)
{

}

void InputOpenFoamParallel::nextVariables(Mesh &mesh)
{
	eslog::start("TIME STEP LOADER: STARTED", "TIME STEP LOADER");
	eslog::info(" == TIME STEP %77s == \n", timesteps[timestep].c_str());

	if (info::mpi::rank == 0) {
		if (utils::exists(variablePath + "/processor" + std::to_string(0) + "/" + timesteps[timestep] + "/uniform/time")) {
			step::time.current = OpenFOAMTime::value(variablePath + "/processor" + std::to_string(0) + "/" + timesteps[timestep] + "/uniform/time");
		}
	}
	Communication::broadcast(&step::time.current, 1, MPI_DOUBLE, 0);

	variablePack.wait(MPITools::singleton->within);

	eslog::checkpointln("TIME STEP LOADER: VARIABLES LOADED");

	for (int d = info::mpi::rank, offset = 0; d < domains; d += info::mpi::size, offset += nvariables.size() + evariables.size()) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			if (variablePack.files[v + offset]->totalSize) {
				OpenFOAMVectorField::load(variablePack.files[v + offset], variables.nodes[v].data);
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			if (variablePack.files[v + nvariables.size() + offset]->totalSize) {
				OpenFOAMVectorField::load(variablePack.files[v + nvariables.size() + offset], variables.elements[v].data);
			}
		}
	}

	if (++timestep < timesteps.size()) {
		variablePack.clear();
		for (int d = info::mpi::rank; d < domains; d += info::mpi::size) {
			for (size_t v = 0; v < nvariables.size(); ++v) {
				variablePack.add(variablePath + "/processor" + std::to_string(d) + "/" + timesteps[timestep] + "/" + variableNames[nvariables[v]]);
			}
			for (size_t v = 0; v < evariables.size(); ++v) {
				variablePack.add(variablePath + "/processor" + std::to_string(d) + "/" + timesteps[timestep] + "/" + variableNames[evariables[v]]);
			}
		}
		variablePack.iread(MPITools::singleton->within);
		eslog::checkpointln("TIME STEP LOADER: IREAD NEXT TIME STEP");
	}

	builder::chunkedValuesNext(variables, mesh);
	eslog::endln("TIME STEP LOADER: FINISHED");
}

void InputOpenFoamParallelDirect::nextVariables(Mesh &mesh)
{

}

void InputOpenFoamSequential::ivariables(const InputConfiguration &configuration)
{

}

void InputOpenFoamParallel::ivariables(const InputConfiguration &configuration)
{
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

	std::vector<std::string> subdirs;
	for (int d = info::mpi::rank; d < info::mpi::size; d += info::mpi::size) { // only once
		utils::listDirectorySubdirectories(configuration.path + "/processor" + std::to_string(d), subdirs);
	}

	std::vector<std::string> steps;
	if (configuration.openfoam.time.size()) {
		std::vector<std::string> values = Parser::split(configuration.openfoam.time, ":");
		if (values.size() == 1) {
			values.push_back(values.front());
		}
		for (size_t s = 0; s < subdirs.size(); ++s) {
			if (StringCompare::caseInsensitiveEq(subdirs[s], values.front())) {
				steps.push_back(subdirs[s]); continue;
			}
			if (StringCompare::caseInsensitiveEq(subdirs[s], values.back())) {
				steps.push_back(subdirs[s]); continue;
			}
			if (StringCompare::caseInsensitive(values.front(), subdirs[s]) && StringCompare::caseInsensitive(subdirs[s], values.back())) {
				steps.push_back(subdirs[s]); continue;
			}
		}
	} else { // get the last subdir
		int max = -1;
		float maxtime = 0;
		for (size_t s = 0; s < subdirs.size(); ++s) {
			char *strend;
			float time = std::strtof(subdirs[s].c_str(), &strend);
			if (max == -1) {
				if (subdirs[s].c_str() != strend) {
					max = s;
					maxtime = time;
				}
			} else {
				if (maxtime < time) {
					max = s;
					maxtime = time;
				}
			}
		}
		if (max != -1) {
			steps.push_back(subdirs[max]);
		}
	}

	if (steps.size() == 0) {
		return;
	}

	std::vector<double> times(steps.size());
	for (size_t t = 0; t < times.size(); ++t) {
		std::stringstream ss(steps[t]);
		ss >> times[t];
	}
	std::vector<int> perm(steps.size());
	std::iota(perm.begin(), perm.end(), 0);
	std::sort(perm.begin(), perm.end(), [&] (int i, int j) { return times[i] < times[j]; });
	for (size_t i = 0; i < perm.size(); ++i) {
		timesteps.push_back(steps[perm[i]]);
	}

	eslog::info(" == TIME STEPS       %70d == \n", timesteps.size());
	eslog::info(" == > FIRST          %70s == \n", timesteps.front().c_str());
	eslog::info(" == > LAST           %70s == \n", timesteps.back().c_str());

	eslog::info(" == VARIABLES %77s == \n", "");
	std::vector<int> skip;
	for (int d = info::mpi::rank; d < info::mpi::size; d += info::mpi::size) { // only once
		if (configuration.openfoam.variables.size()) {
			variableNames = Parser::split(configuration.openfoam.variables, ",");
		} else {
			utils::listDirectoryFiles(configuration.path + "/processor" + std::to_string(d) + "/" + timesteps.back(), variableNames);
		}
		if (utils::exists(configuration.path + "/processor" + std::to_string(d) + "/" + timesteps.back() + "/polyMesh/points")) {
			variableNames.push_back("polyMesh/points");
		}
		skip.resize(variableNames.size(), false);
		vheader.resize(variableNames.size());
		for (size_t i = 0; i < variableNames.size(); ++i) {
			std::ifstream is(configuration.path + "/processor" + std::to_string(d) + "/" + timesteps.back() + "/" + variableNames[i]);
			vheader[i].read(is);
			std::string name;
			switch (vheader[i].foamClass) {
			case FoamFileHeader::Class::pointScalarField:   nvariables.push_back(i); name = "pointScalarField"; break;
			case FoamFileHeader::Class::pointVectorField:   nvariables.push_back(i); name = "pointVectorField"; break;
			case FoamFileHeader::Class::vectorField:        nvariables.push_back(i); name = "vectorField"; break;
			case FoamFileHeader::Class::volScalarField:     evariables.push_back(i); name = "volScalarField"; break;
			case FoamFileHeader::Class::volVectorField:     evariables.push_back(i); name = "volVectorField"; break;
			case FoamFileHeader::Class::surfaceScalarField: svariables.push_back(i); name = "[skipped] surfaceScalarField"; skip[i] = true; break;
			case FoamFileHeader::Class::surfaceVectorField: svariables.push_back(i); name = "[skipped] surfaceVectorField"; skip[i] = true; break;
			default: name = "unknown"; break;
			}
			eslog::info(" == > %s %*s == \n", variableNames[i].c_str(), 84 - variableNames[i].size(), name.c_str());
		}
	}

	variablePath = configuration.path;

	for (int d = info::mpi::rank; d < domains; d += info::mpi::size) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			variablePack.add(variablePath + "/processor" + std::to_string(d) + "/" + timesteps.front() + "/" + variableNames[nvariables[v]]);
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			variablePack.add(variablePath + "/processor" + std::to_string(d) + "/" + timesteps.front() + "/" + variableNames[evariables[v]]);
		}
	}
	variablePack.iread(MPITools::singleton->within);

	for (int d = info::mpi::rank; d < info::mpi::size; d += info::mpi::size) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			if (StringCompare::caseInsensitiveEq("polyMesh/points", variableNames[nvariables[v]])) {
				variables.nodes.push_back({vheader[nvariables[v]].dimension(), "coordinates"});
			} else {
				variables.nodes.push_back({vheader[nvariables[v]].dimension(), variableNames[nvariables[v]]});
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			variables.elements.push_back({vheader[evariables[v]].dimension(), variableNames[evariables[v]]});
		}
	}
}

void InputOpenFoamParallelDirect::ivariables(const InputConfiguration &configuration)
{

}


