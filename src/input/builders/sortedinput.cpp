
#include <config/ecf/ecf.h>
#include "sortedinput.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

SortedInput::SortedInput(MeshBuilder &mesh)
: Input(mesh)
{
	eslog::startln("BUILDER: BUILD SORTED MESH", "BUILDER");

	balance();
	checkERegions();
	eslog::checkpointln("BUILDER: BALANCE DATA");

	fillElements();
	eslog::checkpointln("BUILDER: ELEMENTS FILLED");

	fillCoordinates();
	eslog::checkpointln("BUILDER: NODES FILLED");

	addNodeRegions();
	eslog::checkpointln("BUILDER: NODES REGIONS FILLED");

	addBoundaryRegions();
	eslog::checkpointln("BUILDER: BOUNDARY REGIONS FILLED");

	addElementRegions();
	eslog::endln("BUILDER: ELEMENTS REGIONS FILLED");
}

void SortedInput::checkERegions()
{
//	std::vector<MeshERegion> bregions;
//
//	for (size_t r = 0; r < _meshData.eregions.size(); r++) {
//		if (_meshData.eregions[r].min < _eDistribution.back() && _eDistribution.back() < _meshData.eregions[r].max) {
//			ESINFO(ERROR) << "ESPRESO Workbench parser error: weird element region.";
//		}
//		if (_meshData.eregions[r].min >= _eDistribution.back()) {
//			bregions.push_back(MeshERegion(std::move(_meshData.eregions[r])));
//			_meshData.eregions.erase(_meshData.eregions.begin() + r--);
//		}
//	}
//
//	size_t bsize = 0;
//	std::vector<size_t> rsize = { 0 };
//	for (size_t i = 0; i < _meshData.bregions.size(); i++) {
//		bsize += _meshData.bregions[i].esize.size();
//		rsize.push_back(bsize);
//	}
//
//	std::vector<size_t> fdistribution = Communication::getDistribution(bsize, MPITools::operations().sizeToOffsetsSize_t);
//
//	size_t origBSize = _meshData.bregions.size();
//
//	for (size_t r = 0; r < bregions.size(); r++) {
//		std::vector<size_t> borders;
//		for (int t = 0; t < info::mpi::MPIsize; t++) {
//			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
//			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
//			if (begin != end) {
//				borders.push_back(*begin);
//				borders.push_back(borders.back() + end - begin);
//			}
//		}
//
//		if (!Communication::allGatherUnknownSize(borders)) {
//			ESINFO(ERROR) << "ESPRESO internal error: gather bregion borders.";
//		}
//
//		bool onlyRename = false;
//		for (size_t br = 0; br < origBSize; br++) {
//			if (_meshData.bregions[br].min == borders.front() && _meshData.bregions[br].max == borders.back() - 1) {
//				_meshData.bregions[br].name = bregions[r].name;
//				onlyRename = true;
//				break;
//			}
//		}
//		if (onlyRename) {
//			continue;
//		}
//
//		std::vector<int> tRanks;
//		std::vector<std::vector<esint> > sBuffer, rBuffer;
//
//		for (int t = 0; t < info::mpi::MPIsize; t++) {
//			auto begin = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t] + _eDistribution.back());
//			auto end = std::lower_bound(bregions[r].elements.begin(), bregions[r].elements.end(), fdistribution[t + 1] + _eDistribution.back());
//			if (begin != end) {
//				tRanks.push_back(t);
//				sBuffer.push_back(std::vector<esint>(begin, end));
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: send boundary region indices.";
//		}
//
//		for (size_t i = 1; i < rBuffer.size(); i++) {
//			rBuffer[0].insert(rBuffer[0].end(), rBuffer[i].begin(), rBuffer[i].end());
//		}
//
//		auto cmp = [] (EData &edata, esint id) {
//			return edata.id < id;
//		};
//
//		_meshData.bregions.push_back(MeshBRegion());
//		_meshData.bregions.back().name = bregions[r].name;
//		if (rBuffer.size() && rBuffer.front().size()) {
//			for (size_t nr = 0; nr < origBSize; nr++) {
//				if (_meshData.bregions[nr].esize.size()) {
//					auto begin = std::lower_bound(_meshData.bregions[nr].edata.begin(), _meshData.bregions[nr].edata.end(), rBuffer[0].front(), cmp);
//					auto end = std::lower_bound(_meshData.bregions[nr].edata.begin(), _meshData.bregions[nr].edata.end(), rBuffer[0].back() + 1, cmp);
//					for (size_t i = begin - _meshData.bregions[nr].edata.begin(), nodes = 0; i < end - _meshData.bregions[nr].edata.begin(); nodes += _meshData.bregions[nr].esize[i++]) {
//						_meshData.bregions.back().edata.push_back(_meshData.bregions[nr].edata[i]);
//						_meshData.bregions.back().enodes.insert(_meshData.bregions.back().enodes.end(), _meshData.bregions[nr].enodes.begin() + nodes, _meshData.bregions[nr].enodes.begin() + nodes + _meshData.bregions[nr].esize[i]);
//						_meshData.bregions.back().esize.push_back(_meshData.bregions[nr].esize[i]);
//					}
//				}
//			}
//		}
//	}
}

void SortedInput::fillCoordinates()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	if (info::mpi::size == 1) {
		std::vector<std::vector<Point> > tcoordinates(threads);
		std::vector<std::vector<esint> > nIDs(threads), rData(threads);

		std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, _meshData.coordinates.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			tcoordinates[t].insert(tcoordinates[t].end(), _meshData.coordinates.begin() + cdistribution[t], _meshData.coordinates.begin() + cdistribution[t + 1]);
			nIDs[t].insert(nIDs[t].end(), _meshData.nIDs.begin() + cdistribution[t], _meshData.nIDs.begin() + cdistribution[t + 1]);
		}

		info::mesh->nodes->size = _meshData.coordinates.size();
		info::mesh->nodes->distribution = cdistribution;
		info::mesh->nodes->IDs = new serializededata<esint, esint>(1, nIDs);
		info::mesh->nodes->coordinates = new serializededata<esint, Point>(1, tcoordinates);
		info::mesh->nodes->ranks = new serializededata<esint, int>(1, tarray<int>(threads, _nDistribution.back()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
			std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
		}
		info::mesh->boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
		info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, rData);

		info::mesh->neighborsWithMe.push_back(info::mpi::rank);
		return;
	}

	std::vector<std::vector<esint> > nodes(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tnodes(info::mesh->elements->procNodes->datatarray().begin(t), info::mesh->elements->procNodes->datatarray().end(t));
		utils::sortAndRemoveDuplicates(tnodes);
		nodes[t].swap(tnodes);
	}
	utils::inplaceMerge(nodes);
	utils::removeDuplicates(nodes[0]);

	std::vector<std::vector<esint> > sBuffer;
	std::vector<int> sRanks;
	std::vector<int> ssize(info::mpi::size), rsize(info::mpi::size);

	for (int t = 0; t < info::mpi::size; t++) {
		auto begin = std::lower_bound(nodes[0].begin(), nodes[0].end(), _nDistribution[t]);
		auto end = std::lower_bound(nodes[0].begin(), nodes[0].end(), _nDistribution[t + 1]);
		if (end - begin) {
			sBuffer.push_back(std::vector<esint>(begin, end));
			sRanks.push_back(t);
		}
		ssize[t] = end - begin;
	}

	////////

	std::vector<esint> rrIDs;
	MPI_Alltoall(ssize.data(), 1, MPI_INT, rsize.data(), 1, MPI_INT, info::mpi::comm);

	size_t rrsize = 0;
	for (int t = 0; t < info::mpi::size; t++) {
		rrsize += rsize[t];
	}
	rrIDs.resize(rrsize);

	Communication::allToAllV(nodes[0], rrIDs, ssize, rsize);

	if (!Communication::sendVariousTargets(sBuffer, _rankNodeMap, sRanks, _targetRanks)) {
		eslog::error("ESPRESO internal error: exchange neighbors.\n");
	}
	{
		size_t rnodesize = 0;
		for (size_t t = 0; t < _targetRanks.size(); t++) {
			rnodesize += _rankNodeMap[t].size();
		}
		if (rnodesize != rrIDs.size()) {
			eslog::error("INVALID ALL TO ALL EXCHANGE.\n");
		}
	}

	std::vector<esint> ndistribution = tarray<esint>::distribute(threads, _meshData.coordinates.size());
	std::vector<std::vector<std::vector<esint> > > backedData(threads, std::vector<std::vector<esint> >(_targetRanks.size()));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> ranks, ranksOffset;
		std::vector<std::vector<esint> > tbackedData(_targetRanks.size());
		std::vector<std::vector<esint>::const_iterator> rPointer(_targetRanks.size());

		for (size_t r = 0; r < _targetRanks.size(); r++) {
			rPointer[r] = std::lower_bound(_rankNodeMap[r].begin(), _rankNodeMap[r].end(), _nDistribution[info::mpi::rank] + ndistribution[t]);
		}
		for (esint n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			ranks.clear();
			ranksOffset.clear();
			for (size_t r = 0; r < _targetRanks.size(); r++) {
				if (rPointer[r] != _rankNodeMap[r].end() && *rPointer[r] == _nDistribution[info::mpi::rank] + n) {
					ranksOffset.push_back(r);
					ranks.push_back(_targetRanks[r]);
					++rPointer[r];
				}
			}
			for (size_t r = 0; r < ranks.size(); r++) {
				tbackedData[ranksOffset[r]].push_back(ranksOffset.size());
				tbackedData[ranksOffset[r]].insert(tbackedData[ranksOffset[r]].end(), ranks.begin(), ranks.end());
			}
		}

		backedData[t].swap(tbackedData);
	}

	#pragma omp parallel for
	for (size_t r = 0; r < _targetRanks.size(); r++) {
		for (size_t t = 1; t < threads; t++) {
			backedData[0][r].insert(backedData[0][r].end(), backedData[t][r].begin(), backedData[t][r].end());
		}
	}

	std::vector<std::vector<Point> > backedCoordinates(_targetRanks.size());
	#pragma omp parallel for
	for (size_t r = 0; r < _targetRanks.size(); r++) {
		backedCoordinates[r].resize(_rankNodeMap[r].size());
	}

	for (size_t r = 0; r < _targetRanks.size(); r++) {
		std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, _rankNodeMap[r].size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = rdistribution[t]; n < rdistribution[t + 1]; ++n) {
				backedCoordinates[r][n] = _meshData.coordinates[_rankNodeMap[r][n] - _nDistribution[info::mpi::rank]];
			}
		}
	}

	std::vector<std::vector<esint> > nodeRanks(sRanks.size()), allnodes(threads);
	std::vector<std::vector<Point> > coordinates(sRanks.size());

	if (!Communication::sendVariousTargets(backedData[0], nodeRanks, _targetRanks)) {
		eslog::error("ESPRESO internal error: return node ranks.\n");
	}
	if (!Communication::sendVariousTargets(backedCoordinates, coordinates, _targetRanks)) {
		eslog::error("ESPRESO internal error: return coordinates.\n");
	}

	size_t csize = 0;
	for (size_t i = 0; i < coordinates.size(); i++) {
		csize += coordinates[i].size();
	}

	std::vector<size_t> distribution = tarray<size_t>::distribute(threads, csize);
	std::vector<std::vector<esint> > rankDistribution(sRanks.size());
	std::vector<std::vector<int> > rankData(sRanks.size());

	#pragma omp parallel for
	for (size_t r = 0; r < sRanks.size(); r++) {
		std::vector<esint> trankDistribution;
		std::vector<int> trankData;
		if (r == 0) {
			trankDistribution.push_back(0);
		}

		for (size_t n = 0; n < nodeRanks[r].size(); n += nodeRanks[r][n] + 1) {
			trankData.insert(trankData.end(), nodeRanks[r].begin() + n + 1, nodeRanks[r].begin() + n + 1 + nodeRanks[r][n]);
			trankDistribution.push_back(trankData.size());
		}

		rankDistribution[r].swap(trankDistribution);
		rankData[r].swap(trankData);
	}

	utils::threadDistributionToFullDistribution(rankDistribution);

	for (size_t i = threads; i < sRanks.size(); i++) {
		coordinates[threads - 1].insert(coordinates[threads - 1].end(), coordinates[i].begin(), coordinates[i].end());
		rankData[threads - 1].insert(rankData[threads - 1].end(), rankData[i].begin(), rankData[i].end());
		rankDistribution[threads - 1].insert(rankDistribution[threads - 1].end(), rankDistribution[i].begin(), rankDistribution[i].end());
	}
	for (size_t i = threads; i < sRanks.size(); i++) {
		sBuffer[threads - 1].insert(sBuffer[threads - 1].end(), sBuffer[i].begin(), sBuffer[i].end());
	}
	coordinates.resize(threads);
	sBuffer.resize(threads);
	rankData.resize(threads);
	rankDistribution.resize(threads);

	serializededata<esint, Point>::balance(1, coordinates, &distribution);
	serializededata<esint, esint>::balance(1, sBuffer, &distribution);
	serializededata<esint, int>::balance(rankDistribution, rankData, &distribution);

	info::mesh->nodes->size = distribution.back();
	info::mesh->nodes->distribution = distribution;
	info::mesh->nodes->IDs = new serializededata<esint, esint>(1, sBuffer);
	info::mesh->nodes->coordinates = new serializededata<esint, Point>(1, coordinates);
	info::mesh->nodes->ranks = new serializededata<esint, int>(rankDistribution, rankData);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		allnodes[t].resize(distribution[t + 1] - distribution[t]);
		std::iota(allnodes[t].begin(), allnodes[t].end(), distribution[t]);
		utils::sortAndRemoveDuplicates(rankData[t]);
	}

	info::mesh->boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
	info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, allnodes);

	for (size_t t = 0; t < threads; t++) {
		info::mesh->neighborsWithMe.insert(info::mesh->neighborsWithMe.end(), rankData[t].begin(), rankData[t].end());
	}
	utils::sortAndRemoveDuplicates(info::mesh->neighborsWithMe);

	for (size_t n = 0; n < info::mesh->neighborsWithMe.size(); n++) {
		if (info::mesh->neighborsWithMe[n] != info::mpi::rank) {
			info::mesh->neighbors.push_back(info::mesh->neighborsWithMe[n]);
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = info::mesh->elements->procNodes->begin(t)->begin(); n != info::mesh->elements->procNodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), *n) - info::mesh->nodes->IDs->datatarray().begin();
		}
	}
}

void SortedInput::addNodeRegions()
{
	// assume sorted nodes !!
	size_t threads = info::env::OMP_NUM_THREADS;

	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion) {
		std::vector<std::vector<esint> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < info::mpi::size; t++) {
			auto begin = std::lower_bound(nregion->second.begin(), nregion->second.end(), _nDistribution[t]);
			auto end = std::lower_bound(nregion->second.begin(), nregion->second.end(), _nDistribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<esint>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			eslog::error("ESPRESO internal error: exchange node region.\n");
		}

		sBuffer.clear();
		sBuffer.resize(_targetRanks.size());
		for (size_t r = 1; r < rBuffer.size(); r++) {
			rBuffer[0].insert(rBuffer[0].end(), rBuffer[r].begin(), rBuffer[r].end());
		}

		if (rBuffer.size()) {
			#pragma omp parallel for
			for (size_t t = 0; t < _targetRanks.size(); t++) {
				sBuffer[t].resize(rBuffer[0].size());
				sBuffer[t].resize(std::set_intersection(_rankNodeMap[t].begin(), _rankNodeMap[t].end(), rBuffer[0].begin(), rBuffer[0].end(), sBuffer[t].begin()) - sBuffer[t].begin());
			}
		}

		for (size_t t = 0; t < _targetRanks.size(); t++) {
			if (sBuffer[t].size()) {
				tRanks.push_back(t);
			}
		}
		for (size_t t = 0; t < tRanks.size(); t++) {
			sBuffer[t].swap(sBuffer[tRanks[t]]);
			tRanks[t] = _targetRanks[tRanks[t]];
		}
		sBuffer.resize(tRanks.size());

		rBuffer.clear();
		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			eslog::error("ESPRESO internal error: exchange node region to targets.\n");
		}

		for (size_t t = threads; t < rBuffer.size(); t++) {
			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
		}
		rBuffer.resize(threads);
		serializededata<esint, esint>::balance(1, rBuffer);

		info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(nregion->first));
		info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, rBuffer);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto n = info::mesh->boundaryRegions.back()->nodes->begin(t)->begin(); n != info::mesh->boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
				*n = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), *n) - info::mesh->nodes->IDs->datatarray().begin();
			}
		}
	}
}

void SortedInput::addBoundaryRegions()
{
//	size_t threads = info::env::OMP_NUM_THREADS;
//
//	if (info::mpi::MPIsize == 1) {
//		for (size_t i = 0; i < _meshData.bregions.size(); i++) {
//			std::vector<esint> edist = { 0 };
//			edist.reserve(_meshData.bregions[i].esize.size() + 1);
//			for (size_t e = 0; e < _meshData.bregions[i].esize.size(); e++) {
//				edist.push_back(edist.back() + _meshData.bregions[i].esize[e]);
//			}
//
//			std::vector<std::vector<esint> > tedist(threads), tnodes(threads);
//			std::vector<std::vector<Element*> > epointers(threads);
//			std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _meshData.bregions[i].esize.size());
//
//			tedist.front().push_back(0);
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				for (size_t n = edistribution[t]; n < edistribution[t + 1]; ++n) {
//					tnodes[t].insert(tnodes[t].end(), _meshData.bregions[i].enodes.begin() + edist[n], _meshData.bregions[i].enodes.begin() + edist[n + 1]);
//					epointers[t].push_back(&info::mesh->_eclasses[t][_meshData.bregions[i].edata[n].etype]);
//					tedist[t].push_back(tnodes[t].size());
//				}
//			}
//
//			utils::threadDistributionToFullDistribution(tedist);
//
//			info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(_meshData.bregions[i].name, info::mesh->_eclasses));
//			info::mesh->boundaryRegions.back()->distribution = tarray<esint>::distribute(threads, epointers.front().size());
//			switch (epointers.front().front()->type) {
//			case Element::TYPE::PLANE:
//				info::mesh->boundaryRegions.back()->dimension = 2;
//				break;
//			case Element::TYPE::LINE:
//				info::mesh->boundaryRegions.back()->dimension = 1;
//				break;
//			default:
//				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
//			}
//			info::mesh->boundaryRegions.back()->elements = new serializededata<esint, esint>(tedist, tnodes);
//			info::mesh->boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, epointers);
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				for (auto n = info::mesh->boundaryRegions.back()->elements->begin(t)->begin(); n != info::mesh->boundaryRegions.back()->elements->end(t)->begin(); ++n) {
//					*n = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), *n) - info::mesh->nodes->IDs->datatarray().begin();
//				}
//			}
//		}
//		return;
//	}
//
//	TimeEval timing("BOUNDARY REGIONS");
//	timing.totalTime.startWithBarrier();
//
//	TimeEvent e1("BR LINK NODES AND ELEMENTS"); e1.start();
//
//	if (_meshData.bregions.size()) {
//		info::mesh->preprocessing->linkNodesAndElements();
//	}
//
//	e1.end(); timing.addEvent(e1);
//
//	std::vector<esint> edistribution = info::mesh->elements->gatherElementsProcDistribution();
//
//	for (size_t i = 0; i < _meshData.bregions.size(); i++) {
//
//		TimeEvent e2("BR PREPARE"); e2.start();
//
//		std::vector<esint> edist = { 0 };
//		edist.reserve(_meshData.bregions[i].esize.size() + 1);
//		for (size_t e = 0; e < _meshData.bregions[i].esize.size(); e++) {
//			edist.push_back(edist.back() + _meshData.bregions[i].esize[e]);
//		}
//
//		std::vector<esint> permutation(edist.size() - 1);
//		std::iota(permutation.begin(), permutation.end(), 0);
//		std::sort(permutation.begin(), permutation.end(), [&] (esint e1, esint e2) {
//			return _meshData.bregions[i].enodes[edist[e1]] < _meshData.bregions[i].enodes[edist[e2]];
//		});
//
//		e2.end(); timing.addEvent(e2);
//
//		TimeEvent e3("BR SRANKS"); e3.start();
//
//		std::vector<std::vector<esint> > sBuffer, rBuffer;
//		std::vector<int> sRanks, tRanks;
//
//		for (int t = 0; t < info::mpi::MPIsize; t++) {
//			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t], [&] (esint e, esint n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
//			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t + 1], [&] (esint e, esint n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
//			if (begin != end) {
//				sRanks.push_back(t);
//			}
//		}
//		sBuffer.resize(sRanks.size());
//
//		e3.end(); timing.addEvent(e3);
//
//		TimeEvent e4("BR SBUFFER"); e4.start();
//
//		for (size_t r = 0; r < sRanks.size(); r++) {
//			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[r]], [&] (esint e, esint n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
//			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[r] + 1], [&] (esint e, esint n) { return _meshData.bregions[i].enodes[edist[e]] < n; });
//			std::vector<size_t> sdistribution = tarray<esint>::distribute(threads, end - begin);
//			std::vector<std::vector<esint> > tsBuffer(threads);
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				std::vector<esint> ttsBuffer;
//
//				for (auto e = begin + sdistribution[t]; e != begin + sdistribution[t + 1]; ++e) {
//					ttsBuffer.push_back(_meshData.bregions[i].edata[*e].etype);
//					ttsBuffer.push_back(_meshData.bregions[i].esize[*e]);
//					for (esint n = 0; n < _meshData.bregions[i].esize[*e]; ++n) {
//						ttsBuffer.push_back(_meshData.bregions[i].enodes[edist[*e] + n]);
//					}
//				}
//
//				tsBuffer[t].swap(ttsBuffer);
//			}
//
//			sBuffer[r].push_back(0);
//			for (size_t t = 0; t < threads; t++) {
//				sBuffer[r].push_back(tsBuffer[t].size() + sBuffer[r].back());
//			}
//			for (size_t t = 0; t < threads; t++) {
//				sBuffer[r].insert(sBuffer[r].end(), tsBuffer[t].begin(), tsBuffer[t].end());
//			}
//
//		}
//
//		e4.end(); timing.addEvent(e4);
//
//		int avgneighs = 0, nneighs = sRanks.size();
//		double allavgsize = 0, avgsize = 0;
//		for (size_t j = 0; j < sRanks.size(); j++) {
//			avgsize += sBuffer[j].size();
//		}
//		avgsize /= nneighs;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_SUM, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "BR AVGNEIGHS: " << (double)avgneighs / info::mpi::MPIsize << ", AVGSIZE: " << allavgsize / info::mpi::MPIsize;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MIN, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MIN, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "BR MINNEIGHS: " << avgneighs << ", MINSIZE: " << allavgsize;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "BR MAXNEIGHS: " << avgneighs << ", MAXSIZE: " << allavgsize;
//
//		TimeEvent e5("BR EXCHANGE SBUFFER"); e5.start();
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
//		}
//
//		e5.end(); timing.addEvent(e5);
//
//
//		nneighs = rBuffer.size();
//		avgsize = 0;
//		for (size_t j = 0; j < rBuffer.size(); j++) {
//			avgsize += rBuffer[j].size();
//		}
//		avgsize /= nneighs;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_SUM, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "AVGNEIGHS: " << (double)avgneighs / info::mpi::MPIsize << ", AVGSIZE: " << allavgsize / info::mpi::MPIsize;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MIN, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MIN, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "MINNEIGHS: " << avgneighs << ", MINSIZE: " << allavgsize;
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, info::mpi::MPICommunicator);
//		MPI_Reduce(&avgsize, &allavgsize, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::MPICommunicator);
//
//		ESINFO(PROGRESS1) << "MAXNEIGHS: " << avgneighs << ", MAXSIZE: " << allavgsize;
//
//		nneighs = _targetRanks.size();
//
//		MPI_Reduce(&nneighs, &avgneighs, 1, MPI_INT, MPI_MAX, 0, info::mpi::MPICommunicator);
//		ESINFO(PROGRESS1) << "AVGTARGETS: " << (double)avgneighs / info::mpi::MPIsize;
//
//		TimeEvent e6("BR PROCESS RBUFFER"); e6.start();
//
//		sBuffer.clear();
//		sBuffer.resize(_targetRanks.size());
//
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//			std::vector<std::vector<std::vector<esint> > > tsBuffer(threads, std::vector<std::vector<esint> >(_targetRanks.size()));
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				std::vector<esint> nodes;
//				std::vector<std::vector<esint> > ttsBuffer(_targetRanks.size());
//
//				for (size_t n = rBuffer[r][t] + threads + 1; n < rBuffer[r][t + 1] + threads + 1; n += 2 + rBuffer[r][n + 1]) {
//					nodes.clear();
//					nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//					std::sort(nodes.begin(), nodes.end());
//					auto nbegin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[info::mpi::MPIrank]);
//					auto nend = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[info::mpi::MPIrank + 1]);
//
//					for (size_t tt = 0; tt < _targetRanks.size(); tt++) {
//						auto it = _rankNodeMap[tt].begin();
//						bool found = true;
//						for (auto current = nbegin; found && current != nend; ++current) {
//							it = std::lower_bound(it, _rankNodeMap[tt].end(), *current);
//							found = it != _rankNodeMap[t].end() && *it == *current;
//						}
//						if (found) {
//							ttsBuffer[tt].insert(ttsBuffer[tt].end(), rBuffer[r].begin() + n, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
//						}
//					}
//				}
//
//				tsBuffer[t].swap(ttsBuffer);
//			}
//
//			for (size_t tt = 0; tt < _targetRanks.size(); tt++) {
//				size_t tsize = 0;
//				for (size_t t = 0; t < threads; t++) {
//					tsize += tsBuffer[t][tt].size();
//				}
//				if (tsize) {
//					sBuffer[tt].push_back(0);
//					for (size_t t = 0; t < threads; t++) {
//						sBuffer[tt].push_back(sBuffer[tt].back() + tsBuffer[t][tt].size());
//					}
//				}
//				for (size_t t = 0; t < threads; t++) {
//					sBuffer[tt].insert(sBuffer[tt].end(), tsBuffer[t][tt].begin(), tsBuffer[t][tt].end());
//				}
//			}
//		}
//
//		e6.end(); timing.addEvent(e6);
//
//		TimeEvent e7("BR SEND DATA TO POTENTIAL OWNERS"); e7.start();
//
//		for (size_t t = 0; t < _targetRanks.size(); t++) {
//			if (sBuffer[t].size()) {
//				tRanks.push_back(t);
//			}
//		}
//		for (size_t t = 0; t < tRanks.size(); t++) {
//			sBuffer[t].swap(sBuffer[tRanks[t]]);
//			tRanks[t] = _targetRanks[tRanks[t]];
//		}
//		sBuffer.resize(tRanks.size());
//
//		rBuffer.clear();
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
//		}
//
//		e7.end(); timing.addEvent(e7);
//
//		TimeEvent e8("BR BUILD FACES"); e8.start();
//
//		std::vector<std::vector<esint> > tedist(threads), tnodes(threads);
//		std::vector<std::vector<Element*> > epointers(threads);
//
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				std::vector<esint> ttedist, ttnodes;
//				std::vector<Element*> tepointers;
//				if (t == 0 && r == 0) {
//					ttedist.push_back(0);
//				}
//				esint foffset = 0;
//				if (r && tedist[t].size()) {
//					foffset = tedist[t].back();
//				}
//
//				std::vector<esint> nodes;
//				std::vector<esint> nlinks;
//				int counter;
//				bool found = true;
//				for (esint e = rBuffer[r][t] + threads + 1; e < rBuffer[r][t + 1] + threads + 1; e += 2 + rBuffer[r][e + 1]) {
//					found = true;
//					for (auto n = rBuffer[r].begin() + e + 2; found && n != rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]; ++n) {
//						auto it = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), *n);
//						if (it != info::mesh->nodes->IDs->datatarray().end() && *it == *n) {
//							*n = it - info::mesh->nodes->IDs->datatarray().begin();
//						} else {
//							found = false;
//						}
//					}
//					if (found) {
//						nlinks.clear();
//						for (auto n = rBuffer[r].begin() + e + 2; n != rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]; ++n) {
//							auto links = info::mesh->nodes->elements->cbegin() + *n;
//							nlinks.insert(nlinks.end(), links->begin(), links->end());
//						}
//						std::sort(nlinks.begin(), nlinks.end());
//						counter = 1;
//						for (size_t i = 1; i < nlinks.size(); ++i) {
//							if (nlinks[i - 1] == nlinks[i]) {
//								++counter;
//								if (counter == rBuffer[r][e + 1]) {
//									if (_eDistribution[info::mpi::MPIrank] <= nlinks[i] && nlinks[i] < _eDistribution[info::mpi::MPIrank + 1]) {
//										ttnodes.insert(ttnodes.end(), rBuffer[r].begin() + e + 2, rBuffer[r].begin() + e + 2 + rBuffer[r][e + 1]);
//										ttedist.push_back(ttnodes.size() + foffset);
//										tepointers.push_back(&info::mesh->_eclasses[0][rBuffer[r][e]]);
//									}
//									break;
//								}
//							} else {
//								counter = 1;
//							}
//						}
//					}
//				}
//
//				tedist[t].insert(tedist[t].end(), ttedist.begin(), ttedist.end());
//				tnodes[t].insert(tnodes[t].end(), ttnodes.begin(), ttnodes.end());
//				epointers[t].insert(epointers[t].end(), tepointers.begin(), tepointers.end());
//			}
//		}
//
//		e8.end(); timing.addEvent(e8);
//
//		TimeEvent e10("BR CREATE ARRAYS"); e10.start();
//
//		utils::threadDistributionToFullDistribution(tedist);
//
//		serializededata<esint, esint>::balance(tedist, tnodes);
//		serializededata<esint, Element*>::balance(1, epointers);
//
//		#pragma omp parallel for
//		for (size_t t = 1; t < threads; t++) {
//			for (size_t e = 0; e < epointers[t].size(); e++) {
//				epointers[t][e] = &info::mesh->_eclasses[t][epointers[t][e] - info::mesh->_eclasses[0]];
//			}
//		}
//
//		info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(_meshData.bregions[i].name, info::mesh->_eclasses));
//		info::mesh->boundaryRegions.back()->dimension = 2;
//		if (epointers.front().size()) {
//			switch (epointers.front().front()->type) {
//			case Element::TYPE::PLANE:
//				info::mesh->boundaryRegions.back()->dimension = 2;
//				break;
//			case Element::TYPE::LINE:
//				info::mesh->boundaryRegions.back()->dimension = 1;
//				break;
//			default:
//				ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
//			}
//		}
//		int dim = info::mesh->boundaryRegions.back()->dimension;
//		MPI_Allreduce(&dim, &info::mesh->boundaryRegions.back()->dimension, 1, MPI_INT, MPI_MIN, info::mpi::MPICommunicator);
//
//		info::mesh->boundaryRegions.back()->elements = new serializededata<esint, esint>(tedist, tnodes);
//		info::mesh->boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, epointers);
//		info::mesh->boundaryRegions.back()->distribution = info::mesh->boundaryRegions.back()->epointers->datatarray().distribution();
//
//		e10.end(); timing.addEvent(e10);
//
//		TimeEvent e11("-------"); e11.start();
//		e11.end(); timing.addEvent(e11);
//	}
//
//	timing.totalTime.endWithBarrier();
//	timing.printStatsMPI();
}

void SortedInput::addElementRegions()
{
//	if (info::mpi::MPIsize == 1) {
//		for (size_t i = 0; i < _meshData.eregions.size(); i++) {
//			info::mesh->elementsRegions.push_back(new ElementsRegionStore(_meshData.eregions[i].name));
//
//			std::vector<size_t> distribution = tarray<esint>::distribute(threads, _meshData.eregions[i].elements.size());
//			std::vector<std::vector<esint> > telements(threads);
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				telements[t].insert(telements[t].end(), _meshData.eregions[i].elements.begin() + distribution[t], _meshData.eregions[i].elements.begin() + distribution[t + 1]);
//			}
//			info::mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(1, telements);
//		}
//		return;
//	}
//
//	for (size_t i = 0; i < _meshData.eregions.size(); i++) {
//		std::vector<std::vector<esint> > sBuffer, rBuffer;
//		std::vector<int> sRanks, tRanks;
//
//		for (int t = 0; t < info::mpi::MPIsize; t++) {
//			auto begin = std::lower_bound(_meshData.eregions[i].elements.begin(), _meshData.eregions[i].elements.end(), _eDistribution[t]);
//			auto end = std::lower_bound(_meshData.eregions[i].elements.begin(), _meshData.eregions[i].elements.end(), _eDistribution[t + 1]);
//			if (end - begin) {
//				sBuffer.push_back(std::vector<esint>(begin, end));
//				sRanks.push_back(t);
//			}
//		}
//
//		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
//			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
//		}
//
//		for (size_t t = threads; t < rBuffer.size(); t++) {
//			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
//		}
//		rBuffer.resize(threads);
//		serializededata<esint, esint>::balance(1, rBuffer);
//
//		info::mesh->elementsRegions.push_back(new ElementsRegionStore(_meshData.eregions[i].name));
//		info::mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(1, rBuffer);
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (auto e = info::mesh->elementsRegions.back()->elements->begin(t)->begin(); e != info::mesh->elementsRegions.back()->elements->end(t)->begin(); ++e) {
//				*e -= _eDistribution[info::mpi::MPIrank];
//			}
//		}
//	}
}
