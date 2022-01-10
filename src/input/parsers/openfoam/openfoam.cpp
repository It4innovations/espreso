
#include "openfoam.h"

#include "basis/containers/tarray.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "parser/connectivity.h"
#include "parser/geometry.h"
#include "input/builders/builder.h"

using namespace espreso;

void InputOpenFoam::load(const InputConfiguration &configuration)
{
	eslog::startln("OPENFOAM PARSER: STARTED", "OPENFOAM PARSER");

	InputFilePack geofiles;
	OpenFOAMGeometry geometry(
			geofiles.add(configuration.path + "/constant/polyMesh/points"),
			geofiles.add(configuration.path + "/constant/polyMesh/faces"));

	OpenFOAMConnectivity connectivity(
			mesh.regions->files.add(configuration.path + "/constant/polyMesh/owner"),
			mesh.regions->files.add(configuration.path + "/constant/polyMesh/neighbour"));

	geofiles.prepare();
	mesh.regions->files.setTotalSizes();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READER PREPARED");

	geofiles.read();
	profiler::synccheckpoint("geometry_read");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READ");

	geometry.scan();
	profiler::synccheckpoint("scan");

	// compute better distribution for binary data according to scanned parameters
	connectivity.ownerFile->setDistribution(tarray<size_t>::distribute(info::mpi::size, connectivity.ownerFile->totalSize));
	connectivity.neighborsFile->setDistribution(tarray<size_t>::distribute(info::mpi::size, connectivity.neighborsFile->totalSize));
	mesh.regions->files.iread();
	mesh.regions->files.wait();

	geometry.parse(mesh.nodes, mesh.elements);
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");

	eslog::endln("OPENFOAM PARSER: PARSED");
}

void InputOpenFoam::build(Mesh &mesh)
{
	builder::build(this->mesh, mesh);
}

//#include "openfoam.h"
//
//#include "parser/points.h"
//#include "parser/faces.h"
//#include "parser/boundary.h"
//#include "parser/zones.h"
//#include "parser/sets.h"
//
//#include "input/meshbuilder.h"
//#include "basis/containers/tarray.h"
//#include "wrappers/mpi/communication.h"
//#include "basis/utilities/parser.h"
//#include "basis/utilities/utils.h"
//#include "esinfo/envinfo.h"
//#include "esinfo/mpiinfo.h"
//#include "esinfo/eslog.hpp"
//#include "config/ecf/input/input.h"
//#include "mesh/element.h"
//
//#include <numeric>
//
//using namespace espreso;
//
//std::string OpenFOAMLoader::cellprefix = "cells_";
//
//
//OpenFOAMLoader::OpenFOAMLoader(const InputConfiguration &configuration)
//: _configuration(configuration)
//{
//
//}
//
//void OpenFOAMLoader::load()
//{
//	eslog::start("OPENFOAM: STARTED", "OPENFOAM");
//	eslog::param("database", _configuration.path.c_str());
//	eslog::ln();
//
//	readData();
//	eslog::checkpointln("OPENFOAM: DATA READ");
//
//	parseData();
//	eslog::checkpointln("OPENFOAM: DATA PARSED");
//
//	buildFaces();
//	eslog::checkpointln("OPENFOAM: FACE BUILT");
//
//	collectFaces();
//	eslog::checkpointln("OPENFOAM: FACE COLLECTED");
//
//	buildElements();
//	eslog::checkpointln("OPENFOAM: ELEMENTS BUILT");
//
//	eslog::endln("OPENFOAM: MESH BUILT");
//}
//
//void OpenFOAMLoader::readData()
//{
//	eslog::start("LOADER: READER STARTED", "LOADER");
//	eslog::param("READERS", MPITools::singleton->across.size);
//	eslog::ln();
//
//	eslog::start("LOADER: POINTS", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_points.readBySubset(_configuration.path + "/constant/polyMesh/points");
////	_points.alignlines(0);
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: POINTS READ");
//
//	eslog::start("LOADER: FACES", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_faces.readBySubset(_configuration.path + "/constant/polyMesh/faces");
////	_faces.alignlines(0);
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: FACES READ");
//
//	eslog::start("LOADER: OWNERS", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_owner.readBySubset(_configuration.path + "/constant/polyMesh/owner");
////	_owner.alignlines(0);
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: OWNERS READ");
//
//	eslog::start("LOADER: NEIGHBOURS", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_neighbor.readBySubset(_configuration.path + "/constant/polyMesh/neighbour");
////	_neighbor.alignlines(0);
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: NEIGHBOURS READ");
//
//	eslog::start("LOADER: BOUNDARY", "LOADER"); eslog::param("READERS", MPITools::singleton->across.size); eslog::ln();
////	_boundary.readBySingleton(_configuration.path + "/constant/polyMesh/boundary");
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: BOUNDARY READ");
//
//	eslog::start("LOADER: POINTS ZONES", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_pointZones.readBySubset(_configuration.path + "/constant/polyMesh/pointZones");
////	if (_pointZones.totalsize) {
////		_pointZones.alignlines(0);
////	}
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: POINTS ZONES READ (IF ANY)");
//
//	eslog::start("LOADER: FACE ZONES", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_faceZones.readBySubset(_configuration.path + "/constant/polyMesh/faceZones");
////	if (_faceZones.totalsize) {
////		_faceZones.alignlines(0);
////	}
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: FACE ZONES READ (IF ANY)");
//
//	eslog::start("LOADER: CELL ZONES", "LOADER"); eslog::param("READERS", MPITools::subset->across.size); eslog::ln();
////	_cellZones.readBySubset(_configuration.path + "/constant/polyMesh/cellZones");
////	if (_cellZones.totalsize) {
////		_cellZones.alignlines(0);
////	}
//	eslog::endln("LOADER: FINISHED");
//	eslog::checkpointln("OPENFOAM: FACE CELL READ (IF ANY)");
//
//	if (info::mpi::rank == 0) {
//		OpenFOAMSets::inspect(_configuration.path + "/constant/polyMesh/sets/*", _sets);
//		size_t size = _sets.size();
//		Communication::broadcast(&size, sizeof(size_t), MPI_BYTE, 0);
//		Communication::broadcast(_sets.data(), _sets.size() * sizeof(OpenFOAMSet), MPI_BYTE, 0);
//	} else {
//		size_t size;
//		Communication::broadcast(&size, sizeof(size_t), MPI_BYTE, 0);
//		_sets.resize(size);
//		Communication::broadcast(_sets.data(), _sets.size() * sizeof(OpenFOAMSet), MPI_BYTE, 0);
//	}
//
//	eslog::endln("OPENFOAM: PARSED");
//}
//
//void OpenFOAMLoader::parseData()
//{
//	if (!OpenFOAMPoints(_points.begin, _points.end).readData(nIDs, coordinates, 1)) {
//		eslog::error("OpenFOAM loader: cannot parse points.\n");
//	}
//	if (!OpenFOAMFaces(_faces.begin, _faces.end).readFaces(*this)) {
//		eslog::error("OpenFOAM loader: cannot parse faces.\n");
//	}
//	if (!OpenFOAMFaces(_owner.begin, _owner.end).readParents(this->owner)) {
//		eslog::error("OpenFOAM loader: cannot parse owners.\n");
//	}
//	if (!OpenFOAMFaces(_neighbor.begin, _neighbor.end).readParents(this->neighbor)) {
//		eslog::error("OpenFOAM loader: cannot parse neighbors.\n");
//	}
//
//	if (!OpenFOAMBoundary(_boundary.begin, _boundary.end).readData(*this)) {
//		eslog::error("OpenFOAM loader: cannot parse boundary.\n");
//	}
//
//	if (_pointZones.distribution.back() != 0 && !OpenFOAMZones(_pointZones).readPoints(*this)) {
//		eslog::error("OpenFOAM loader: cannot parse pointZones.\n");
//	}
//
//	if (_faceZones.distribution.back() != 0 && !OpenFOAMZones(_faceZones).readFaces(*this)) {
//		eslog::error("OpenFOAM loader: cannot parse faceZones.\n");
//	}
//
//	if (_cellZones.distribution.back() != 0 && !OpenFOAMZones(_cellZones).readCells(*this)) {
//		eslog::error("OpenFOAM loader: cannot parse cellZones.\n");
//	}
//}
//
//void OpenFOAMLoader::buildFaces()
//{
//	// 1. Find MAX element ID in order to be able correctly set face IDs in continuous interval
//
//	esint maxID = 0;
//	if (owner.size()) {
//		maxID = *std::max_element(owner.begin(), owner.end());
//	}
//	if (neighbor.size()) {
//		maxID = std::max(*std::max_element(neighbor.begin(), neighbor.end()), maxID);
//	}
//
//	Communication::allReduce(&maxID, &nelements, 1, MPITools::getType<esint>().mpitype, MPI_MAX);
//	nelements += 1;
//
//	_edist = tarray<esint>::distribute(info::mpi::size, nelements);
//
//	std::vector<esint> fDistribution = Communication::getDistribution<esint>(fsize.size());
//
//	// 2. Add sets that are not in any zone or boundary
//
//	for (size_t i = 0; i < _sets.size(); i++) {
//		std::string name = std::string(_sets[i].name);
//		switch (_sets[i].type) {
//		case OpenFOAMSet::SetType::CELL_SET:
//			for (auto ereg = eregions.begin(); ereg != eregions.end(); ++ereg) {
//				if (StringCompare::caseInsensitivePreffix(ereg->first, cellprefix + name)) {
//					_sets.erase(_sets.begin() + i--);
//					break;
//				}
//			}
//			break;
//		case OpenFOAMSet::SetType::FACE_SET:
//			for (auto ereg = eregions.begin(); ereg != eregions.end(); ++ereg) {
//				if (!StringCompare::caseInsensitivePreffix(ereg->first, cellprefix + name)) {
//					_sets.erase(_sets.begin() + i--);
//					break;
//				}
//			}
//			break;
//		case OpenFOAMSet::SetType::POINT_SET:
//			for (auto nreg = nregions.begin(); nreg != nregions.end(); ++nreg) {
//				if (StringCompare::caseInsensitivePreffix(nreg->first, name)) {
//					_sets.erase(_sets.begin() + i--);
//					break;
//				}
//			}
//			break;
//		}
//	}
//
//	eslog::start("LOADER: SETS", "LOADER");
//	eslog::param("READERS", MPITools::subset->across.size);
//	eslog::ln();
//	for (size_t i = 0; i < _sets.size(); i++) {
//		std::string name = std::string(_sets[i].name);
//		switch (_sets[i].type) {
//		case OpenFOAMSet::SetType::CELL_SET: {
//			InputFile file;
////			file.readBySubset(_configuration.path + "/constant/polyMesh/sets/" + name);
////			file.alignlines(0);
//			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], eregions[cellprefix + name])) {
//				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
//			}
//		} break;
//		case OpenFOAMSet::SetType::FACE_SET: {
//			InputFile file;
////			file.readBySubset(_configuration.path + "/constant/polyMesh/sets/" + name);
////			file.alignlines(0);
//			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], eregions[name])) {
//				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
//			}
//		} break;
//		case OpenFOAMSet::SetType::POINT_SET: {
//			InputFile file;
////			file.readBySubset(_configuration.path + "/constant/polyMesh/sets/" + name);
////			file.alignlines(0);
//			if (!OpenFOAMSets(file.begin, file.end).readData(_sets[i], nregions[name])) {
//				eslog::error("OpenFOAM loader: cannot parse set '%s'.\n", name.c_str());
//			}
//		} break;
//		}
//	}
//	eslog::endln("LOADER: FINISHED");
//
//	// 3. Exchange region data to processes that hold given faces
//
//	std::vector<esint> sBuffer, rBuffer;
//
//	std::vector<size_t> rpointer(eregions.size());
//	for (int r = 0; r < info::mpi::size; r++) {
//		size_t prevsize = sBuffer.size();
//		sBuffer.push_back(0); // total size
//		sBuffer.push_back(r); // target
//
//		size_t rindex = 0;
//		for (auto it = eregions.begin(); it != eregions.end(); ++it, ++rindex) {
//			if (!StringCompare::caseInsensitivePreffix(cellprefix, it->first)) {
//				size_t prevrsize = sBuffer.size();
//				sBuffer.push_back(0); // region size
//
//				for ( ; rpointer[rindex] < it->second.size() && it->second[rpointer[rindex]] < fDistribution[r + 1]; ++rpointer[rindex]) {
//					sBuffer.push_back(it->second[rpointer[rindex]]);
//				}
//				sBuffer[prevrsize] = sBuffer.size() - prevrsize - 1;
//			}
//		}
//
//		sBuffer[prevsize] = sBuffer.size() - prevsize;
//	}
//
//	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
//		eslog::internalFailure("distribute faces indices in regions.\n");
//	}
//
//	for (auto it = eregions.begin(); it != eregions.end(); ++it) {
//		if (!StringCompare::caseInsensitivePreffix(cellprefix, it->first)) {
//			it->second.clear();
//		}
//	}
//
//	std::vector<esint> usedfaces;
//	size_t offset = 0;
//	for (int r = 0; r < info::mpi::size; r++) {
//		++offset; // skip total size
//		++offset; // skip target
//
//		size_t rindex = 0;
//		for (auto it = eregions.begin(); it != eregions.end(); ++it, ++rindex) {
//			if (!StringCompare::caseInsensitivePreffix(cellprefix, it->first)) {
//				size_t rsize = rBuffer[offset++];
//
//				for (size_t i = 0; i < rsize; ++i) {
//					usedfaces.push_back(rBuffer[offset]);
//					it->second.push_back(rBuffer[offset++]); // index faces after elements
//				}
//			}
//			std::sort(it->second.begin(), it->second.end());
//		}
//	}
//
//	// 4. Add used faces into elements
//
//	utils::sortAndRemoveDuplicates(usedfaces);
//
//	size_t foffset = usedfaces.size();
//	Communication::exscan(foffset);
//	foffset += nelements;
//
//	_fdist.reserve(fsize.size() + 1);
//	_fdist.push_back(0);
//	for (size_t f = 0; f < fsize.size(); f++) {
//		_fdist.push_back(_fdist.back() + fsize[f]);
//	}
//
//	for (size_t i = 0; i < usedfaces.size(); i++) {
//		esint findex = usedfaces[i] - fDistribution[info::mpi::rank];
//
//		esize.push_back(fsize[findex]);
//		enodes.insert(enodes.end(), fnodes.begin() + _fdist[findex], fnodes.begin() + _fdist[findex + 1]);
//		if (fsize[findex] == 3) {
//			etype.push_back((int)Element::CODE::TRIANGLE3);
//		}
//		if (fsize[findex] == 4) {
//			etype.push_back((int)Element::CODE::SQUARE4);
//		}
//	}
//
//	eIDs.resize(esize.size(), 0);
//	std::iota(eIDs.begin(), eIDs.end(), foffset);
//	body.resize(esize.size(), 0);
//	material.resize(esize.size(), 0);
//
//	for (auto it = eregions.begin(); it != eregions.end(); ++it) {
//		if (!StringCompare::caseInsensitivePreffix(cellprefix, it->first)) {
//			auto id = usedfaces.begin();
//			for (size_t i = 0; i < it->second.size(); i++) {
//				while (*id < it->second[i]) { ++id; }
//				it->second[i] = foffset + id - usedfaces.begin();
//			}
//		}
//	}
//}
//
//void OpenFOAMLoader::collectFaces()
//{
//	std::vector<size_t> ownersDist = Communication::getDistribution(owner.size());
//	std::vector<size_t> neighborsDist = Communication::getDistribution(neighbor.size());
//	std::vector<size_t> target = Communication::getDistribution(fsize.size());
//
//	if (!Communication::balance(owner, ownersDist, target)) {
//		eslog::internalFailure("balance faces owners.\n");
//	}
//
//	for (size_t i = 0; i < target.size(); i++) {
//		if (target[i] > neighborsDist.back()) {
//			target[i] = neighborsDist.back();
//		}
//	}
//
//	if (!Communication::balance(neighbor, neighborsDist, target)) {
//		eslog::internalFailure("balance faces neighbors.\n");
//	}
//
//	size_t firstID = Communication::getDistribution(fsize.size())[info::mpi::rank];
//
//	auto sortIDs = [&] (std::vector<esint> &permutation, const std::vector<esint> &data) {
//		permutation.resize(data.size());
//		std::iota(permutation.begin(), permutation.end(), 0);
//
//		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
//			if (data[i] == data[j]) {
//				return i < j;
//			}
//			return data[i] < data[j];
//		});
//	};
//
//	std::vector<esint> oPermutation, nPermutation;
//	sortIDs(oPermutation, owner);
//	sortIDs(nPermutation, neighbor);
//
//	std::vector<esint> sBuffer, rBuffer;
//	// ID, size, owner / -1 * neighbor, nodes
//	sBuffer.reserve(4 * info::mpi::size + 3 * (owner.size() + neighbor.size()) + fnodes.size());
//
//	size_t prevsize;
//	auto obegin = oPermutation.begin();
//	auto nbegin = nPermutation.begin();
//	for (int r = 0; r < info::mpi::size; r++) {
//		prevsize = sBuffer.size();
//		sBuffer.push_back(0); // total size
//		sBuffer.push_back(r); // target
//		sBuffer.push_back(0); // number of faces
//
//		auto o = obegin;
//		for ( ; o != oPermutation.end() && owner[*o] < _edist[r + 1]; ++o) {
//			sBuffer.push_back(firstID + *o);
//			sBuffer.push_back(owner[*o]);
//			sBuffer.push_back(fsize[*o]);
//			sBuffer.insert(sBuffer.end(), fnodes.begin() + _fdist[*o], fnodes.begin() + _fdist[*o + 1]);
//		}
//		auto n = nbegin;
//		for ( ; n != nPermutation.end() && neighbor[*n] < _edist[r + 1]; ++n) {
//			sBuffer.push_back(firstID + *n);
//			sBuffer.push_back(-neighbor[*n] - 1);
//			sBuffer.push_back(fsize[*n]);
//			sBuffer.insert(sBuffer.end(), fnodes.begin() + _fdist[*n], fnodes.begin() + _fdist[*n + 1]);
//		}
//		sBuffer[prevsize + 2] = (o - obegin) + (n - nbegin);
//		obegin = o;
//		nbegin = n;
//
//		sBuffer[prevsize] = sBuffer.size() - prevsize;
//	}
//
//	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
//		eslog::internalFailure("distribute permuted elements.\n");
//	}
//
//	fIDs.clear();
//	fsize.clear();
//	fnodes.clear();
//	owner.clear();
//	neighbor.clear();
//
//	std::vector<esint> fIDs, fsize, fnodes, owners;
//
//	size_t offset = 0;
//	for (int r = 0; r < info::mpi::size; r++) {
//		++offset;
//		size_t size = rBuffer[++offset];
//		++offset;
//
//		for (size_t f = 0; f < size; ++f) {
//			fIDs.push_back(rBuffer[offset++]);
//			owners.push_back(rBuffer[offset++]);
//			fsize.push_back(rBuffer[offset++]);
//			fnodes.insert(fnodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + fsize.back());
//			offset += fsize.back();
//		}
//	}
//
//	std::vector<esint> fpermutation(fIDs.size());
//	std::iota(fpermutation.begin(), fpermutation.end(), 0);
//	std::sort(fpermutation.begin(), fpermutation.end(), [&] (esint i, esint j) {
//		return fIDs[i] < fIDs[j];
//	});
//
//	_fdist.resize(1);
//	_fdist.reserve(fsize.size() + 1);
//	for (size_t e = 0; e < fsize.size(); e++) {
//		_fdist.push_back(_fdist.back() + fsize[e]);
//	}
//
//	if (fpermutation.size()) {
//		size_t rest = 0;
//		for (auto i = fpermutation.begin(); i != fpermutation.end() - 1; ++i) {
//			if (fIDs[*i] == fIDs[*(i + 1)]) {
//				fsize.push_back(fsize[*i]);
//				if (owners[*i] >= 0) {
//					owner.push_back(owners[*i]);
//					neighbor.push_back(-owners[*(i + 1)] - 1);
//					fnodes.insert(fnodes.end(), fnodes.begin() + _fdist[*i], fnodes.begin() + _fdist[*i + 1]);
//				} else {
//					owner.push_back(owners[*(i + 1)]);
//					neighbor.push_back(-owners[*i] - 1);
//					fnodes.insert(fnodes.end(), fnodes.begin() + _fdist[*(i + 1)], fnodes.begin() + _fdist[*(i + 1) + 1]);
//				}
//				++i;
//			} else {
//				fpermutation[rest++] = *i;
//			}
//		}
//		if (fIDs[*(fpermutation.end() - 1)] != fIDs[*(fpermutation.end() - 2)]) {
//			fpermutation[rest++] = fpermutation.back();
//		}
//
//		for (auto i = fpermutation.begin(); i != fpermutation.begin() + rest; ++i) {
//			fsize.push_back(fsize[*i]);
//			if (owners[*i] >= 0) {
//				owner.push_back(owners[*i]);
//				fnodes.insert(fnodes.end(), fnodes.begin() + _fdist[*i], fnodes.begin() + _fdist[*i + 1]);
//			} else {
//				owner.push_back(-owners[*i] - 1);
//				fnodes.insert(fnodes.end(), fnodes.rbegin() + _fdist.back() - _fdist[*i + 1], fnodes.rbegin() + _fdist.back() - _fdist[*i]);
//			}
//		}
//	}
//}
//
//void OpenFOAMLoader::buildElements()
//{
//	size_t threads = info::env::OMP_NUM_THREADS;
//
//	auto sortIDs = [&] (std::vector<esint> &permutation, const std::vector<esint> &data) {
//		permutation.resize(data.size());
//		std::iota(permutation.begin(), permutation.end(), 0);
//
//		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
//			if (data[i] == data[j]) {
//				return i < j;
//			}
//			return data[i] < data[j];
//		});
//	};
//
//	std::vector<esint> owner, neighbor;
//	sortIDs(owner, owner);
//	sortIDs(neighbor, neighbor);
//
//	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, _edist[info::mpi::rank + 1] - _edist[info::mpi::rank]);
//
//	std::vector<std::vector<esint> > tesize(threads), tenodes(threads);
//	std::vector<std::vector<int> > tetype(threads);
//
//	_fdist.clear();
//	_fdist.reserve(tdistribution.back() + 1);
//	_fdist.push_back(0);
//	for (size_t f = 0; f < fsize.size(); f++) {
//		_fdist.push_back(_fdist.back() + fsize[f]);
//	}
//
//	auto getThreadBegin = [] (const std::vector<esint> &data, const std::vector<esint> &perm, esint eindex) {
//		return std::lower_bound(perm.begin(), perm.end(), eindex, [&] (esint i, esint eindex) {
//			return data[i] < eindex;
//		}) - perm.begin();
//	};
//
//	auto addFaces = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t &triangles, size_t &squares, size_t &index, esint element) {
//		while (index < perm.size() && data[perm[index]] == element) {
//			switch (fsize[perm[index++]]) {
//			case 3: ++triangles; break;
//			case 4:   ++squares; break;
//			}
//		}
//	};
//
//	auto getFace = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t index, esint element, esint n1, esint n2) {
//		while (index < perm.size() && data[perm[index]] == element) {
//			for (esint f = 0; f < fsize[perm[index]]; f++) {
//				if (n1 == fnodes[_fdist[perm[index]] + f] && n2 == fnodes[_fdist[perm[index]] + (f + 1) % fsize[perm[index]]]) {
//					return std::pair<size_t, esint>(index, f);
//				}
//			}
//			++index;
//		}
//		return std::pair<size_t, esint>(index, -1);
//	};
//
//	auto getUnknown = [&] (esint *kbegin, esint *kend, esint *ubegin, esint *uend) {
//		for (auto i = ubegin, j = kbegin; i != uend; ++i) {
//			for (j = kbegin; j != kend; ++j) {
//				if (*i == *j) {
//					break;
//				}
//			}
//			if (j == kend) {
//				return *i;
//			}
//		}
//		return (esint)-1;
//	};
//
//	auto findElementWithSize = [&] (const std::vector<esint> &perm, size_t &index, size_t &max, int size) {
//		while (index < max) { // there is at least one owner
//			if (fsize[perm[index]] == size) {
//				break;
//			} else {
//				++index;
//			}
//		}
//	};
//
//	size_t eoffset = _edist[info::mpi::rank];
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		std::vector<esint> tsize, tnodes;
//		std::vector<int> ttype;
//
//		size_t oindex = getThreadBegin(owner, owner, tdistribution[t] + eoffset);
//		size_t nindex = getThreadBegin(neighbor, neighbor, tdistribution[t] + eoffset);
//
//		for (size_t e = tdistribution[t] + eoffset; e < tdistribution[t + 1] + eoffset; e++) {
//			size_t obegin = oindex, nbegin = nindex;
//			std::pair<size_t, esint> index;
//			size_t triangles = 0, squares = 0;
//			addFaces(owner, owner, triangles, squares, oindex, e);
//			addFaces(neighbor, neighbor, triangles, squares, nindex, e);
//
//			if (squares == 6 && triangles == 0) {
//				size_t ebegin = tnodes.size();
//				ttype.push_back((int)Element::CODE::HEXA8);
//				tsize.push_back(8);
//				tnodes.insert(tnodes.end(), 8, -1);
//				if (obegin < oindex) { // there is at least one owner
//					tnodes[ebegin + 0] = fnodes[_fdist[owner[obegin]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[owner[obegin]] + 1];
//					tnodes[ebegin + 4] = fnodes[_fdist[owner[obegin]] + 3];
//					tnodes[ebegin + 5] = fnodes[_fdist[owner[obegin]] + 2];
//					++obegin;
//				} else {
//					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nbegin]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nbegin]] + 3];
//					tnodes[ebegin + 4] = fnodes[_fdist[neighbor[nbegin]] + 2];
//					tnodes[ebegin + 5] = fnodes[_fdist[neighbor[nbegin]] + 1];
//					++nbegin;
//				}
//
//				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
//				if (index.first < oindex) {
//					tnodes[ebegin + 2] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
//					tnodes[ebegin + 3] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
//				} else {
//					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
//					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
//					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
//				}
//				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 4], tnodes[ebegin + 5]);
//				if (index.first < oindex) {
//					tnodes[ebegin + 6] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
//					tnodes[ebegin + 7] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
//				} else {
//					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin + 5], tnodes[ebegin + 4]);
//					tnodes[ebegin + 6] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
//					tnodes[ebegin + 7] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
//				}
//				continue;
//			}
//
//			if (squares == 0 && triangles == 4) {
//				size_t ebegin = tnodes.size();
//				ttype.push_back((int)Element::CODE::TETRA4);
//				tsize.push_back(4);
//				tnodes.insert(tnodes.end(), 4, -1);
//				if (obegin < oindex) {
//					tnodes[ebegin + 0] = fnodes[_fdist[owner[obegin]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[owner[obegin]] + 2];
//					tnodes[ebegin + 2] = fnodes[_fdist[owner[obegin]] + 1];
//					++obegin;
//				} else {
//					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nbegin]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nbegin]] + 1];
//					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[nbegin]] + 2];
//					++nbegin;
//				}
//
//				if (obegin < oindex) {
//					tnodes[ebegin + 3] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
//							fnodes.data() + _fdist[owner[obegin]], fnodes.data() + _fdist[owner[obegin] + 1]);
//				} else {
//					tnodes[ebegin + 3] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
//							fnodes.data() + _fdist[neighbor[nbegin]], fnodes.data() + _fdist[neighbor[nbegin] + 1]);
//				}
//				continue;
//			}
//
//			if (squares == 3 && triangles == 2) {
//				size_t ebegin = tnodes.size();
//				ttype.push_back((int)Element::CODE::PRISMA6);
//				tsize.push_back(6);
//				tnodes.insert(tnodes.end(), 6, -1);
//				size_t otria = obegin, ntria = nbegin;
//				findElementWithSize(owner, otria, oindex, 3);
//				findElementWithSize(neighbor, ntria, nindex, 3);
//
//				if (otria < oindex) {
//					tnodes[ebegin + 0] = fnodes[_fdist[owner[otria]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[owner[otria]] + 1];
//					tnodes[ebegin + 2] = fnodes[_fdist[owner[otria]] + 2];
//					findElementWithSize(owner, ++otria, oindex, 3);
//				} else {
//					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[ntria]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[ntria]] + 2];
//					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[ntria]] + 1];
//					findElementWithSize(neighbor, ++ntria, nindex, 3);
//				}
//
//				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
//				if (fsize[owner[index.first]] == 3) {
//					index = getFace(owner, owner, ++obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
//				}
//				if (index.first < oindex) {
//					tnodes[ebegin + 3] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
//					tnodes[ebegin + 4] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
//				} else {
//					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
//					if (fsize[neighbor[index.first]] == 3) {
//						index = getFace(neighbor, neighbor, ++nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
//					}
//					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
//					tnodes[ebegin + 4] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
//				}
//
//				if (otria < oindex) {
//					tnodes[ebegin + 5] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
//							fnodes.data() + _fdist[owner[otria]], fnodes.data() + _fdist[owner[otria] + 1]);
//				} else {
//					tnodes[ebegin + 5] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
//							fnodes.data() + _fdist[neighbor[ntria]], fnodes.data() + _fdist[neighbor[ntria] + 1]);
//				}
//				continue;
//			}
//
//			if (squares == 1 && triangles == 4) {
//				size_t ebegin = tnodes.size();
//				ttype.push_back((int)Element::CODE::PYRAMID5);
//				tsize.push_back(5);
//				tnodes.insert(tnodes.end(), 5, -1);
//				size_t osquare = obegin, nsquare = nbegin, otria = obegin, ntria = nbegin;
//				findElementWithSize(owner, osquare, oindex, 4);
//				findElementWithSize(neighbor, nsquare, nindex, 4);
//				findElementWithSize(owner, otria, oindex, 3);
//				findElementWithSize(neighbor, ntria, nindex, 3);
//
//				if (osquare < oindex) {
//					tnodes[ebegin + 0] = fnodes[_fdist[owner[osquare]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[owner[osquare]] + 1];
//					tnodes[ebegin + 2] = fnodes[_fdist[owner[osquare]] + 2];
//					tnodes[ebegin + 3] = fnodes[_fdist[owner[osquare]] + 3];
//				} else {
//					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nsquare]] + 0];
//					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nsquare]] + 3];
//					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[nsquare]] + 2];
//					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[nsquare]] + 1];
//				}
//
//				if (otria < oindex) {
//					tnodes[ebegin + 4] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
//							fnodes.data() + _fdist[owner[otria]], fnodes.data() + _fdist[owner[otria] + 1]);
//				} else {
//					tnodes[ebegin + 4] = getUnknown(
//							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
//							fnodes.data() + _fdist[neighbor[ntria]], fnodes.data() + _fdist[neighbor[ntria] + 1]);
//				}
//				continue;
//			}
//			eslog::error("OpenFOAM parser: an unknown element type with '%ld' triangles and '%ld' squares [ID='%ld'].\n", triangles, squares, e);
//		}
//
//		tesize[t].swap(tsize);
//		tenodes[t].swap(tnodes);
//		tetype[t].swap(ttype);
//	}
//
//	for (size_t t = 0; t < threads; t++) {
//		esize.insert(esize.end(), tesize[t].begin(), tesize[t].end());
//		enodes.insert(enodes.end(), tenodes[t].begin(), tenodes[t].end());
//		etype.insert(etype.end(), tetype[t].begin(), tetype[t].end());
//	}
//	body.resize(esize.size(), 0);
//	material.resize(esize.size(), 0);
//	size_t fsize = eIDs.size();
//	eIDs.resize(esize.size());
//	std::iota(eIDs.begin() + fsize, eIDs.end(), eoffset);
//}
//
//
