
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

	this->faces.eblocks.blocks.push_back(DatabaseOffset{0, 0, numberOfCells(owner.header)});

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
	std::vector<esint> nsize(domains / info::mpi::size + 1), esize(domains / info::mpi::size + 1);
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		std::string dict = configuration.path + "/processor" + std::to_string(d) + "/constant/polyMesh/";
		OpenFOAMFaceList::load(dict + "faces", faces.ftype, faces.fnodes, nodes.coordinates.size());
		nsize[i] = OpenFOAMVectorField::load(dict + "points", nodes.coordinates);

		esize[i] = numberOfCells(OpenFOAMLabelList::load(dict + "owner", faces.owner, faces.elements));
		OpenFOAMLabelList::load(dict + "neighbour", faces.neighbor, faces.elements);
		faces.neighbor.resize(faces.owner.size(), -1);
		faces.elements += esize[i];
	}

	ndistribution = Communication::getDistribution((esint)nodes.coordinates.size());
	edistribution = Communication::getDistribution(faces.elements);
	esint nsum = 0, esum = 0;
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		nodes.blocks.push_back(DatabaseOffset{ndistribution[info::mpi::rank] + nsum, nsum, nsize[i]});
		nsum += nsize[i];
		faces.eblocks.blocks.push_back(DatabaseOffset{edistribution[info::mpi::rank] + esum, esum, esize[i]});
		esum += esize[i];
	}
	nblocks = nodes.blocks;
	eblocks = faces.eblocks.blocks;

	// from local to global addressing
	for (size_t f = 0, foffset = 0; f < faces.ftype.size(); foffset += Element::encode(faces.ftype[f++]).nodes) {
		PolyElement poly(faces.ftype[f], faces.fnodes.data() + foffset);
		for (int n = 0; n < Element::encode(faces.ftype[f]).nodes; ++n) {
			if (poly.isNode(n)) {
				faces.fnodes[n + foffset] += ndistribution[info::mpi::rank];
			}
		}
	}
	for (size_t n = 0; n < faces.owner.size(); ++n) {
		faces.owner[n] += edistribution[info::mpi::rank];
	}
	for (size_t n = 0; n < faces.neighbor.size(); ++n) {
		if (faces.neighbor[n] != -1) {
			faces.neighbor[n] += edistribution[info::mpi::rank];
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

int InputOpenFoam::variables()
{
	return loader->timesteps.size();
}

double InputOpenFoam::nextVariables(Mesh &mesh)
{
	return loader->nextVariables(mesh);
}

double InputOpenFoamSequential::nextVariables(Mesh &mesh)
{
	return 0;
}

double InputOpenFoamParallel::nextVariables(Mesh &mesh)
{
	eslog::start("OPENFOAM VARIABLES LOADER: STARTED", "VARIABLES LOADER");
	eslog::param("TIME-STEP", timestep);
	eslog::ln();

	eslog::info(" == TIME STEP %77s == \n", timesteps[timestep].c_str());

	std::vector<esint> sBuffer, rBuffer, rmap(info::mpi::size);

	if (timestep == 0) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			if (StringCompare::caseInsensitiveEq("points", vheader[nvariables[v]].object)) {
				vdata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, "TRANSLATION")); // temporary store for coordinates
			} else {
				switch (vheader[nvariables[v]].dimension()) {
				case 1: vdata.push_back(mesh.nodes->appendData(1, NamedData::DataType::SCALAR, vheader[nvariables[v]].object)); break;
				case 3: vdata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, vheader[nvariables[v]].object)); break;
				default: break;
				}
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			switch (vheader[evariables[v]].dimension()) {
			case 1: vdata.push_back(mesh.elements->appendData(1, NamedData::DataType::SCALAR, vheader[evariables[v]].object)); break;
			case 3: vdata.push_back(mesh.elements->appendData(3, NamedData::DataType::VECTOR, vheader[evariables[v]].object)); break;
			default: break;
			}
		}

		nperm.resize(mesh.nodes->size);
		eperm.resize(mesh.elements->distribution.process.size);
		std::iota(nperm.begin(), nperm.end(), 0);
		std::sort(nperm.begin(), nperm.end(), [&] (esint i, esint j) { return (mesh.nodes->inputOffset->begin() + i)->front() < (mesh.nodes->inputOffset->begin() + j)->front(); });
		std::iota(eperm.begin(), eperm.end(), 0);
		std::sort(eperm.begin(), eperm.end(), [&] (esint i, esint j) { return (mesh.elements->inputOffset->begin() + i)->front() < (mesh.elements->inputOffset->begin() + j)->front(); });

		sBuffer.reserve(info::mpi::size * 5 + nperm.size() + eperm.size());
		auto nit = nperm.begin(), eit = eperm.begin();
		for (int t = 0; t < info::mpi::size; ++t) {
			auto nbegin = nit, ebegin = eit;
			while (nit != nperm.end() && (mesh.nodes->inputOffset->begin() + *nit)->front() < ndistribution[t + 1]) { ++nit; }
			while (eit != eperm.end() && (mesh.elements->inputOffset->begin() + *eit)->front() < edistribution[t + 1]) { ++eit; }

			sBuffer.push_back(5 + (nit - nbegin) + (eit - ebegin)); // size
			sBuffer.push_back(t); // target
			sBuffer.push_back(info::mpi::rank); // source
			sBuffer.push_back(nit - nbegin); // nodes
			sBuffer.push_back(eit - ebegin); // elements
			for (auto it = nbegin; it != nit; ++it) {
				sBuffer.push_back((mesh.nodes->inputOffset->begin() + *it)->front());
			}
			for (auto it = ebegin; it != eit; ++it) {
				sBuffer.push_back((mesh.elements->inputOffset->begin() + *it)->front());
			}
		}

		if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, variableRequests)) {
			eslog::error("cannot exchange output offsets.\n");
		}

		variableRequestsMap.resize(info::mpi::size);
		for (size_t i = 0; i < variableRequests.size(); i += variableRequests[i]) {
			variableRequestsMap[variableRequests[i + 2]] = i;
		}
		sBuffer.clear();
		eslog::checkpointln("VARIABLES LOADER: PERMUTATION EXCHANGED");
	}

	variablePack.wait(MPITools::singleton->within);

	eslog::checkpointln("VARIABLES LOADER: VARIABLES LOADED");

	std::vector<std::vector<esfloat> > data(nvariables.size() + evariables.size());
	for (int d = info::mpi::rank, offset = 0; d < domains; d += info::mpi::size, offset += nvariables.size() + evariables.size()) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			if (variablePack.files[v + offset]->totalSize) {
				OpenFOAMVectorField::load(variablePack.files[v + offset], data[v]);
			} else {
				data[v].resize(vheader[nvariables[v]].dimension() * mesh.nodes->size);
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			if (variablePack.files[v + nvariables.size() + offset]->totalSize) {
				OpenFOAMVectorField::load(variablePack.files[v + nvariables.size() + offset], data[v + nvariables.size()]);
			} else {
				data[v + nvariables.size()].resize(vheader[evariables[v]].dimension() * mesh.elements->distribution.process.size);
			}
		}
	}

	esint datasize = 0, nvsize = 0, evsize = 0;
	for (size_t v = 0; v < nvariables.size(); ++v) {
		nvsize += vheader[nvariables[v]].dimension();
	}
	for (size_t v = 0; v < evariables.size(); ++v) {
		evsize += vheader[evariables[v]].dimension();
	}
	for (size_t i = 0; i < variableRequests.size(); i += variableRequests[i]) {
		datasize += utils::reinterpret_size<esint, esfloat>(nvsize * variableRequests[i + 3]);
		datasize += utils::reinterpret_size<esint, esfloat>(evsize * variableRequests[i + 4]);
	}
	sBuffer.reserve(info::mpi::size * 5 + datasize);
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = variableRequestsMap[r];
		esint nsize = variableRequests[i + 3];
		esint esize = variableRequests[i + 4];
		sBuffer.push_back(5 + utils::reinterpret_size<esint, esfloat>(nvsize * nsize + evsize * esize));
		sBuffer.push_back(r);
		sBuffer.push_back(info::mpi::rank);
		sBuffer.push_back(nsize);
		sBuffer.push_back(esize);
		size_t size = sBuffer.size(); // it is useless if we called correctly sBuffer.reserve(...)
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, esfloat>(nvsize * nsize + evsize * esize));
		esfloat *c = reinterpret_cast<esfloat*>(sBuffer.data() + size);
		for (size_t v = 0; v < nvariables.size(); ++v) {
			for (esint n = 0; n < nsize; ++n) {
				for (int d = 0; d < vheader[nvariables[v]].dimension(); ++d) {
					*c++ = data[v][vheader[nvariables[v]].dimension() * (variableRequests[i + 5 + n] - ndistribution[info::mpi::rank]) + d];
				}
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			for (esint e = 0; e < esize; ++e) {
				for (int d = 0; d < vheader[evariables[v]].dimension(); ++d) {
					*c++ = data[v + nvariables.size()][vheader[evariables[v]].dimension() * (variableRequests[i + 5 + e + nsize] - edistribution[info::mpi::rank]) + d];
				}
			}
		}
	}

	rBuffer.clear();
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("cannot exchange output data.\n");
	}
	sBuffer.clear();
	eslog::checkpointln("VARIABLES LOADER: VARIABLES EXCHANGED");

	for (size_t i = 0; i < rBuffer.size(); i += rBuffer[i]) {
		rmap[rBuffer[i + 2]] = i;
	}

	auto nit = nperm.begin(), eit = eperm.begin();
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = rmap[r];
		esint nsize = rBuffer[i + 3];
		esint esize = rBuffer[i + 4];
		esfloat *c = reinterpret_cast<esfloat*>(rBuffer.data() + i + 5);
		for (size_t v = 0; v < nvariables.size(); ++v) {
			for (esint n = 0; n < nsize; ++n) {
				for (int d = 0; d < vheader[nvariables[v]].dimension(); ++d) {
					vdata[v]->data[*(nit + n) * vheader[nvariables[v]].dimension() + d] = *c++;
				}
			}
		}
		nit += nsize;
		for (size_t v = 0; v < evariables.size(); ++v) {
			for (esint e = 0; e < esize; ++e) {
				for (int d = 0; d < vheader[evariables[v]].dimension(); ++d) {
					vdata[v + nvariables.size()]->data[*(eit + e) * vheader[evariables[v]].dimension() + d] = *c++;
				}
			}
		}
		eit += esize;
	}
	eslog::checkpointln("VARIABLES LOADER: VARIABLES ASSIGNED");

	for (size_t v = 0; v < nvariables.size(); ++v) {
		if (StringCompare::caseInsensitiveEq("points", vheader[nvariables[v]].object)) {
			if (variablePack.files[v]->totalSize) {
				if (mesh.nodes->originCoordinates == nullptr) {
					mesh.nodes->originCoordinates = mesh.nodes->coordinates;
					mesh.nodes->coordinates = new serializededata<esint, Point>(*mesh.nodes->originCoordinates);
				}
				for (size_t n = 0; n < mesh.nodes->coordinates->datatarray().size(); ++n) {
					mesh.nodes->coordinates->datatarray()[n].x = vdata[v]->data[3 * n + 0];
					mesh.nodes->coordinates->datatarray()[n].y = vdata[v]->data[3 * n + 1];
					mesh.nodes->coordinates->datatarray()[n].z = vdata[v]->data[3 * n + 2];
					vdata[v]->data[3 * n + 0] = mesh.nodes->coordinates->datatarray()[n].x - mesh.nodes->originCoordinates->datatarray()[n].x;
					vdata[v]->data[3 * n + 1] = mesh.nodes->coordinates->datatarray()[n].y - mesh.nodes->originCoordinates->datatarray()[n].y;
					vdata[v]->data[3 * n + 2] = mesh.nodes->coordinates->datatarray()[n].z - mesh.nodes->originCoordinates->datatarray()[n].z;
				}
				mesh.updateMesh();
			} else {
				std::fill(vdata[v]->data.begin(), vdata[v]->data.end(), 0);
			}
			eslog::checkpointln("VARIABLES LOADER: COORDINATES UPDATED");
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
		eslog::checkpointln("VARIABLES LOADER: IREAD NEXT TIME STEP");
	}

	if (step::step.substep == 0) {
		info::mesh->output->updateMonitors();
	}
	info::mesh->output->updateSolution();
	eslog::checkpointln("VARIABLES LOADER: SOLUTION UPDATED");

	double time = 0;
	for (int d = info::mpi::rank; d < info::mpi::size; d += info::mpi::size) { // only once
		if (utils::exists(variablePath + "/processor" + std::to_string(d) + "/" + timesteps[timestep - 1] + "/uniform/time")) {
			time = OpenFOAMTime::value(variablePath + "/processor" + std::to_string(d) + "/" + timesteps[timestep - 1] + "/uniform/time");
		}
	}

	eslog::endln("VARIABLES LOADER: TIME SET");

	return time;
}

double InputOpenFoamParallelDirect::nextVariables(Mesh &mesh)
{
	return 0;
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
		double min, max;
		std::stringstream smin(values[0]); smin >> min;
		max = min;
		if (values.size() == 2) {
			std::stringstream smax(values[1]); smax >> max;
		}
		for (size_t s = 0; s < subdirs.size(); ++s) {
			char *strend;
			double time = std::strtod(subdirs[s].c_str(), &strend);
			if (subdirs[s].c_str() != strend) {
				if (min - 1e-9 <= time && time <= max + 1e-9) {
					steps.push_back(subdirs[s]);
				}
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
}

void InputOpenFoamParallelDirect::ivariables(const InputConfiguration &configuration)
{

}


