
#include "openfoam.h"

#include "parser/faceList.h"
#include "parser/labelList.h"
#include "parser/vectorField.h"

#include "basis/containers/tarray.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/sysutils.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/builders/builder.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "wrappers/mpi/communication.h"

#include <fstream>
#include <numeric>
#include <algorithm>

using namespace espreso;

InputOpenFoam::InputOpenFoam()
: loader(nullptr)
{
	ParallelFoamFile::init();
}

InputOpenFoam::~InputOpenFoam()
{
	if (loader) delete loader;
	ParallelFoamFile::finish();
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
	ParallelFoamFile::synchronize({ &points, &faces, &owner, &neighbour });
	profiler::synccheckpoint("scan");

	this->faces.elements.blocks.push_back(DatabaseOffset{0, 0, numberOfCells(owner.header)});

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

	ndistribution = Communication::getDistribution((esint)nodes.coordinates.size());
	edistribution = Communication::getDistribution(localCells);
	esint nsum = 0, esum = 0;
	for (int d = info::mpi::rank, i = 0; d < domains; d += info::mpi::size, ++i) {
		nodes.blocks.push_back(DatabaseOffset{ndistribution[info::mpi::rank] + nsum, nsum, nsize[i]});
		nsum += nsize[i];
		faces.elements.blocks.push_back(DatabaseOffset{edistribution[info::mpi::rank] + esum, esum, esize[i]});
		esum += esize[i];
	}
	nblocks = nodes.blocks;
	eblocks = faces.elements.blocks;

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

void InputOpenFoam::variables(Mesh &mesh)
{
	loader->variables(mesh);
}

void InputOpenFoamSequential::variables(Mesh &mesh)
{

}

void InputOpenFoamParallel::variables(Mesh &mesh)
{
	eslog::startln("OPENFOAM VARIABLES LOADER: STARTED", "VARIABLES LOADER");
	ivector<esint> nperm(mesh.nodes->size), eperm(mesh.elements->distribution.process.size);
	std::iota(nperm.begin(), nperm.end(), 0);
	std::sort(nperm.begin(), nperm.end(), [&] (esint i, esint j) { return (mesh.nodes->outputOffset->begin() + i)->front() < (mesh.nodes->outputOffset->begin() + j)->front(); });
	std::iota(eperm.begin(), eperm.end(), 0);
	std::sort(eperm.begin(), eperm.end(), [&] (esint i, esint j) { return (mesh.elements->outputOffset->begin() + i)->front() < (mesh.elements->outputOffset->begin() + j)->front(); });

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(info::mpi::size * 5 + nperm.size() + eperm.size());

	auto nit = nperm.begin(), eit = eperm.begin();
	for (int t = 0; t < info::mpi::size; ++t) {
		auto nbegin = nit, ebegin = eit;
		while (nit != nperm.end() && (mesh.nodes->outputOffset->begin() + *nit)->front() < ndistribution[t + 1]) { ++nit; }
		while (eit != eperm.end() && (mesh.elements->outputOffset->begin() + *eit)->front() < edistribution[t + 1]) { ++eit; }

		sBuffer.push_back(5 + (nit - nbegin) + (eit - ebegin)); // size
		sBuffer.push_back(t); // target
		sBuffer.push_back(info::mpi::rank); // source
		sBuffer.push_back(nit - nbegin); // nodes
		sBuffer.push_back(eit - ebegin); // elements
		for (auto it = nbegin; it != nit; ++it) {
			sBuffer.push_back((mesh.nodes->outputOffset->begin() + *it)->front());
		}
		for (auto it = ebegin; it != eit; ++it) {
			sBuffer.push_back((mesh.elements->outputOffset->begin() + *it)->front());
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("cannot exchange output offsets.\n");
	}

	esint datasize = 0, nvsize = 0, evsize = 0;
	ivector<esint> rmap(info::mpi::size);
	for (size_t v = 0; v < nvariables.size(); ++v) {
		nvsize += vheader[nvariables[v]].dimension();
	}
	for (size_t v = 0; v < evariables.size(); ++v) {
		evsize += vheader[evariables[v]].dimension();
	}
	for (size_t i = 0; i < rBuffer.size(); i += rBuffer[i]) {
		rmap[rBuffer[i + 2]] = i;
		datasize += utils::reinterpret_size<esint, esfloat>(nvsize * rBuffer[i + 3]);
		datasize += utils::reinterpret_size<esint, esfloat>(evsize * rBuffer[i + 4]);
	}
	sBuffer.clear();
	sBuffer.reserve(info::mpi::size * 5 + datasize);

	eslog::checkpointln("VARIABLES LOADER: PERMUTATION EXCHANGED");

	variablePack.wait(MPITools::singleton->across);

	eslog::checkpointln("VARIABLES LOADER: VARIABLES LOADED");

	std::vector<std::vector<esfloat> > data(nvariables.size() + evariables.size());
	for (int d = info::mpi::rank, offset = 0; d < domains; d += info::mpi::size, offset += nvariables.size() + evariables.size()) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			OpenFOAMVectorField::load(variablePack.files[v + offset], data[v]);
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			OpenFOAMVectorField::load(variablePack.files[v + nvariables.size() + offset], data[v + nvariables.size()]);
		}
	}

	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = rmap[r];
		esint nsize = rBuffer[i + 3];
		esint esize = rBuffer[i + 4];
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
					*c++ = data[v][vheader[nvariables[v]].dimension() * (rBuffer[i + 5 + n] - ndistribution[info::mpi::rank]) + d];
				}
			}
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			for (esint e = 0; e < esize; ++e) {
				for (int d = 0; d < vheader[evariables[v]].dimension(); ++d) {
					*c++ = data[v + nvariables.size()][vheader[evariables[v]].dimension() * (rBuffer[i + 5 + e + nsize] - edistribution[info::mpi::rank]) + d];
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

	std::vector<NamedData*> vdata;
	for (size_t v = 0; v < nvariables.size(); ++v) {
		switch (vheader[nvariables[v]].dimension()) {
		case 1: vdata.push_back(mesh.nodes->appendData(1, NamedData::DataType::SCALAR, vheader[nvariables[v]].object)); break;
		case 3: vdata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, vheader[nvariables[v]].object)); break;
		default: break;
		}
	}
	for (size_t v = 0; v < evariables.size(); ++v) {
		switch (vheader[evariables[v]].dimension()) {
		case 1: vdata.push_back(mesh.elements->appendData(1, NamedData::DataType::SCALAR, vheader[evariables[v]].object)); break;
		case 3: vdata.push_back(mesh.elements->appendData(3, NamedData::DataType::VECTOR, vheader[evariables[v]].object)); break;
		default: break;
		}
	}

	for (size_t i = 0; i < rBuffer.size(); i += rBuffer[i]) {
		rmap[rBuffer[i + 2]] = i;
	}

	nit = nperm.begin(), eit = eperm.begin();
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

	eslog::endln("VARIABLES LOADER: VARIABLES ASSIGNED");
}

void InputOpenFoamSequential::ivariables(const InputConfiguration &configuration)
{

}

void InputOpenFoamParallel::ivariables(const InputConfiguration &configuration)
{
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

	std::string subdir;
	std::vector<std::string> files;
	std::vector<int> skip;
	for (int d = info::mpi::rank; d < info::mpi::size; d += info::mpi::size) { // only once
		std::vector<std::string> subdirs;
		utils::listDirectorySubdirectories(configuration.path + "/processor" + std::to_string(d), subdirs);

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
		if (max == -1) {
			break;
		}
		subdir = subdirs[max];

		eslog::info(" == VARIABLES %77s == \n", "");
		std::string fullSubdir = (configuration.path + "/processor*/" + subdir + "/");
		eslog::info(" == %s%*s == \n", fullSubdir.c_str(), 87 - fullSubdir.size(), "");
		utils::listDirectoryFiles(configuration.path + "/processor" + std::to_string(d) + "/" + subdir, files);
		skip.resize(files.size(), false);
		vheader.resize(files.size());
		for (size_t i = 0; i < files.size(); ++i) {
			std::ifstream is(configuration.path + "/processor" + std::to_string(d) + "/" + subdir + "/" + files[i]);
			vheader[i].read(is);
			std::string name;
			switch (vheader[i].foamClass) {
			case FoamFileHeader::Class::pointScalarField:   nvariables.push_back(i); name = "pointScalarField"; break;
			case FoamFileHeader::Class::pointVectorField:   nvariables.push_back(i); name = "pointVectorField"; break;
			case FoamFileHeader::Class::volScalarField:     evariables.push_back(i); name = "volScalarField"; break;
			case FoamFileHeader::Class::volVectorField:     evariables.push_back(i); name = "volVectorField"; break;
			case FoamFileHeader::Class::surfaceScalarField: svariables.push_back(i); name = "[skipped] surfaceScalarField"; skip[i] = true; break;
			default: name = "unknown"; break;
			}
			eslog::info(" == %*s%s%*s == \n", fullSubdir.size(), "", files[i].c_str(), 87 - files[i].size() - fullSubdir.size(), name.c_str());
		}
	}

	for (int d = info::mpi::rank; d < domains; d += info::mpi::size) {
		for (size_t v = 0; v < nvariables.size(); ++v) {
			variablePack.add(configuration.path + "/processor" + std::to_string(d) + "/" + subdir + "/" + files[nvariables[v]]);
		}
		for (size_t v = 0; v < evariables.size(); ++v) {
			variablePack.add(configuration.path + "/processor" + std::to_string(d) + "/" + subdir + "/" + files[evariables[v]]);
		}
	}
	variablePack.iread(MPITools::singleton->across);
}
