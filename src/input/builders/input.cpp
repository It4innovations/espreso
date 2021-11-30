
#include "input.h"
#include "sequentialinput.h"
#include "scatteredinput.h"

#include "basis/containers/serializededata.h"
#include "basis/structures/kdtree.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <algorithm>
#include <numeric>
#include <set>

#include "basis/utilities/print.h"

using namespace espreso;

void Input::clip()
{
	if (info::ecf->input.clipping_box.apply) {
		std::vector<esint> nids, rdist = { 0 };
		std::vector<int> rdata;
		std::vector<_Point<esfloat> > coordinates;
		std::vector<esint> nstatus(_meshData.coordinates.size(), -1);
		const ClippingBox &box = info::ecf->input.clipping_box;
		for (size_t n = 0; n < _meshData.coordinates.size(); ++n) {
			if (!(
					_meshData.coordinates[n].x < box.min[0] || box.max[0] < _meshData.coordinates[n].x ||
					_meshData.coordinates[n].y < box.min[1] || box.max[1] < _meshData.coordinates[n].y ||
					_meshData.coordinates[n].z < box.min[2] || box.max[2] < _meshData.coordinates[n].z)) {

				nids.push_back(_meshData.nIDs[n]);
				coordinates.push_back(_meshData.coordinates[n]);
			} else {
				removed.insert(_meshData.nIDs[n]);
			}
		}

		_meshData.nIDs.swap(nids);
		_meshData.coordinates.swap(coordinates);
		_meshData.nregions.clear();

		std::vector<esint> enodes, esize, eids;
		std::vector<int> etype, mat, body;
		for (size_t e = 0, offset = 0, noffset; e < _meshData.esize.size(); ++e) {
			bool clip = Mesh::edata[_meshData.etype[e]].type != Element::TYPE::VOLUME;
			noffset = enodes.size();
			for (esint n = 0; n < _meshData.esize[e]; ++n) {
				if (removed.count(_meshData.enodes[n + offset])) {
					clip = true;
					break;
				}
				enodes.push_back(_meshData.enodes[n + offset]);
			}
			if (clip) {
				enodes.resize(noffset);
			} else {
				eids.push_back(_meshData.eIDs[e]);
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				mat.push_back(_meshData.material[e]);
				body.push_back(_meshData.body[e]);
			}
			offset += _meshData.esize[e];
		}

		_meshData.eIDs.swap(eids);
		_meshData.esize.swap(esize);
		_meshData.etype.swap(etype);
		_meshData.enodes.swap(enodes);
		_meshData.material.swap(mat);
		_meshData.body.swap(body);
		_meshData.eregions.clear();

		profiler::synccheckpoint("geometry_clipped");
		eslog::checkpointln("BUILDER: GEOMETRY CLIPPED");
	}
}

void Input::serialize()
{
	if (_meshData.nIDs.size() != _meshData.coordinates.size()) {
		eslog::internalFailure("the size of MeshData::nIDs != MeshData::coordinates.\n");
	}


	_meshData.nodes.set(_meshData.nIDs.size(), 0, 0);

	if (std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end())) {
//		#pragma omp parallel for
		for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
			_meshData.nodes[n].offset = _meshData.nIDs[n];
			_meshData.nodes[n].coordinate = _meshData.coordinates[n];
		}
	} else {
		eslog::internalFailure("Implement a loader for non-sorted nodes.\n");
	}

	std::vector<esint>().swap(_meshData.nIDs);
	std::vector<_Point<esfloat> >().swap(_meshData.coordinates);

	if (_meshData.eIDs.size() != _meshData.esize.size()) {
		eslog::internalFailure("the size of MeshData::eIDs != MeshData::esize.\n");
	}
	if (_meshData.eIDs.size() != _meshData.etype.size()) {
		eslog::internalFailure("the size of MeshData::eIDs != MeshData::etype.\n");
	}

	_meshData.elements.set(_meshData.esize, 0, 0);

	if (std::is_sorted(_meshData.eIDs.begin(), _meshData.eIDs.end())) {
//		#pragma omp parallel for
		for (size_t e = 0, edist = 0; e < _meshData.eIDs.size(); edist += _meshData.esize[e++]) {
			_meshData.elements[e].offset = _meshData.eIDs[e];
			_meshData.elements[e].type = _meshData.etype[e];
			_meshData.elements[e].nodes = _meshData.esize[e];
			for (int n = 0; n < _meshData.elements[e].nodes; ++ n) {
				*_meshData.elements[e].node = _meshData.enodes[edist + n];
			}
		}
	} else {
		eslog::internalFailure("Implement a loader for non-sorted elements.\n");
	}

	std::vector<esint>().swap(_meshData.eIDs);
	std::vector<int>().swap(_meshData.etype);
	std::vector<int>().swap(_meshData.esize);
	std::vector<esint>().swap(_meshData.enodes);
}

void Input::connect()
{
	size_t nregions = _meshData.nregions.size();
	size_t eregions = _meshData.eregions.size();
	size_t nvalues = _meshData.ndata.size();
	size_t evalues = _meshData.edata.size();

	std::vector<esint> rintervals, voffsets, vsum(nregions * nvalues + eregions * evalues);
	rintervals.reserve(2 * (nregions + eregions));
	voffsets.reserve(nregions * nvalues + eregions * evalues);
	for (auto reg = _meshData.nregions.begin(); reg != _meshData.nregions.end(); ++reg) {
		if (reg->second.size() == 2) {
			rintervals.push_back(reg->second[0]);
			rintervals.push_back(reg->second[1]);
		} else if (reg->second.size() != 0) {
			eslog::internalFailure("Implement a loader for non-intervaled regions.\n");
		}
		for (auto val = _meshData.ndata.begin(); val != _meshData.ndata.end(); ++val) {
			esint values = 0;
			auto it = val->second.find(reg->first);
			if (it != val->second.end()) {
				values = it->second.size();
			}
			voffsets.push_back(values);
		}
	}
	for (auto reg = _meshData.eregions.begin(); reg != _meshData.eregions.end(); ++reg) {
		if (reg->second.size() == 2) {
			rintervals.push_back(reg->second[0]);
			rintervals.push_back(reg->second[1]);
		} else if (reg->second.size() != 0) {
			eslog::internalFailure("Implement a loader for non-intervaled regions.\n");
		}
		for (auto val = _meshData.edata.begin(); val != _meshData.edata.end(); ++val) {
			esint values = 0;
			auto it = val->second.find(reg->first);
			if (it != val->second.end()) {
				values = it->second.size();
			}
			voffsets.push_back(values);
		}
	}

	Communication::allReduce(rintervals.data(), NULL, rintervals.size(), MPITools::getType<esint>().mpitype, MPI_MAX);
	Communication::exscan(vsum, voffsets);

	std::cout << rintervals;
	std::cout << voffsets;

//	size_t nvalues = _meshData..size();
}

void Input::balance()
{
	profiler::syncstart("balance");
	if (false) {
		auto isSorted = [] (const std::vector<esint> &ids) {
			int sorted, allSorted;

			sorted = std::is_sorted(ids.begin(), ids.end());
			Communication::allReduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN);

			if (!allSorted) {
				return allSorted;
			}

			esint prev = -1, my = ids.size() ? ids.back() : -1;
			if (info::mpi::rank % 2 == 0) {
				if (info::mpi::rank + 1 < info::mpi::size) {
					MPI_Send(&my, sizeof(esint), MPI_BYTE, info::mpi::rank + 1, 0, info::mpi::comm);
				}
				if (info::mpi::rank) {
					MPI_Recv(&prev, sizeof(esint), MPI_BYTE, info::mpi::rank - 1, 0, info::mpi::comm, MPI_STATUS_IGNORE);
				}
			} else {
				MPI_Recv(&prev, sizeof(esint), MPI_BYTE, info::mpi::rank - 1, 0, info::mpi::comm, MPI_STATUS_IGNORE);
				if (info::mpi::rank + 1 < info::mpi::size) {
					MPI_Send(&my, sizeof(esint), MPI_BYTE, info::mpi::rank + 1, 0, info::mpi::comm);
				}
			}

			if (ids.size()) {
				sorted = prev < ids.front();
			}
			Communication::allReduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN);
			return allSorted;
		};

		if (isSorted(_meshData.nIDs)) {
			balanceNodes();
		} else {
			balancePermutedNodes();
			sortNodes();
		}

		if (isSorted(_meshData.eIDs)) {
			balanceElements();
		} else {
			balancePermutedElements();
			std::vector<esint> permutation(_meshData.esize.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });
			sortElements(permutation);
		}
		return;
	}

	std::vector<esint> edist;
	edist.reserve(_meshData.esize.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<esint> npermutation(_meshData.nIDs.size()), epermutation(_meshData.eIDs.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::iota(epermutation.begin(), epermutation.end(), 0);
	if (!std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end())) {
		std::sort(npermutation.begin(), npermutation.end(), [&] (esint i, esint j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });
	}
	if (!std::is_sorted(_meshData.eIDs.begin(), _meshData.eIDs.end())) {
		std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });
	}
	profiler::synccheckpoint("sort");
	profiler::syncparam("nodes", _meshData.nIDs.size());
	profiler::syncparam("elements", _meshData.eIDs.size());
	profiler::syncparam("enodes", _meshData.enodes.size());

//	Communication::computeSplitters(_meshData.nIDs, npermutation, _nDistribution);
//	Communication::computeSplitters(_meshData.eIDs, epermutation, _eDistribution);
	profiler::synccheckpoint("compute_splitters");

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(
			info::mpi::size * 4 +
			_meshData.nIDs.size() * (1 + sizeof(Point) / sizeof(esint)) +
			_meshData.eIDs.size() * 5 + _meshData.esize.size() + // ID, size, body, material, type
			info::mpi::size // enodes size
			);
	auto nit = npermutation.begin();
	auto eit = epermutation.begin();
	for (int r = 0; r < info::mpi::size; ++r) {
		auto nbegin = nit;
		auto ebegin = eit;
		while (nit != npermutation.end() && _meshData.nIDs[*nit] < _nDistribution[r + 1]) { ++nit; }
		while (eit != epermutation.end() && _meshData.eIDs[*eit] < _eDistribution[r + 1]) { ++eit; }

		size_t tsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(nit - nbegin); // nsize
		sBuffer.push_back(eit - ebegin); // esize

		for (auto n = nbegin; n != nit; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
		}
		for (auto n = nbegin; n != nit; ++n) {
			sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<esint*>(_meshData.coordinates.data() + *n + 1));
		}

		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.push_back(_meshData.eIDs[*e]);
		}
		esint sum = 0;
		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.push_back(_meshData.esize[*e]);
			sum += _meshData.esize[*e];
		}
		sBuffer.push_back(sum);
		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + edist[*e], _meshData.enodes.begin() + edist[*e + 1]);
		}
		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.push_back(_meshData.etype[*e]);
		}
		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.push_back(_meshData.material[*e]);
		}
		for (auto e = ebegin; e != eit; ++e) {
			sBuffer.push_back(_meshData.body[*e]);
		}
		sBuffer[tsize] = sBuffer.size() - tsize;
	}

	profiler::synccheckpoint("sbuffer");
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("distribute permuted nodes.\n");
	}
	profiler::synccheckpoint("exchange");

	_meshData.nIDs.clear();
	_meshData.coordinates.clear();

	_meshData.eIDs.clear();
	_meshData.esize.clear();
	_meshData.enodes.clear();
	_meshData.etype.clear();
	_meshData.material.clear();
	_meshData.body.clear();

	size_t offset = 0;
	for (int r = 0; r < info::mpi::size; ++r) {
		++offset; // total size
		++offset; // target
		esint nsize = rBuffer[offset++]; // nsize
		esint esize = rBuffer[offset++]; // esize

		_meshData.nIDs.insert(_meshData.nIDs.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + nsize);
		offset += nsize;
		_meshData.coordinates.insert(_meshData.coordinates.end(), reinterpret_cast<Point*>(rBuffer.data() + offset), reinterpret_cast<Point*>(rBuffer.data() + offset + nsize * sizeof(Point) / sizeof(esint)));
		offset += nsize * sizeof(Point) / sizeof(esint);

		_meshData.eIDs.insert(_meshData.eIDs.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + esize);
		offset += esize;
		_meshData.esize.insert(_meshData.esize.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + esize);
		offset += esize;
		esint sum = rBuffer[offset++];
		_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + sum);
		offset += sum;
		_meshData.etype.insert(_meshData.etype.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + esize);
		offset += esize;
		_meshData.material.insert(_meshData.material.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + esize);
		offset += esize;
		_meshData.body.insert(_meshData.body.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + esize);
		offset += esize;
	}
	profiler::synccheckpoint("rbuffer");

	sortNodes();
	epermutation.resize(_meshData.eIDs.size());
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });
	sortElements(epermutation);
	profiler::synccheckpoint("sort");
	profiler::syncparam("nodes", _meshData.nIDs.size());
	profiler::syncparam("elements", _meshData.eIDs.size());
	profiler::syncparam("enodes", _meshData.enodes.size());
	profiler::syncend("balance");
}

void Input::balanceNodes()
{
	std::vector<esint> cCurrent = Communication::getDistribution<esint>(_meshData.nIDs.size());
	_nDistribution = tarray<esint>::distribute(info::mpi::size, cCurrent.back());

	if (!Communication::balance(_meshData.nIDs, cCurrent, _nDistribution)) {
		eslog::internalFailure("balance node IDs.\n");
	}
	if (!Communication::balance(_meshData.coordinates, cCurrent, _nDistribution)) {
		eslog::internalFailure("balance coordinates.\n");
	}

	esint back = 0, max;
	if (_meshData.nIDs.size()) {
		back = _meshData.nIDs.back();
	}
	Communication::allReduce(&back, &max, 1, MPITools::getType<esint>().mpitype, MPI_MAX);
	if (_meshData.nIDs.size()) {
		back = _meshData.nIDs.back();
	} else {
		back = max;
	}
	Communication::allGather(&back, _nDistribution.data() + 1, sizeof(esint), MPI_BYTE);
	for (size_t i = 1; i < _nDistribution.size(); i++) {
		++_nDistribution[i];
	}
}

void Input::balancePermutedNodes()
{
	std::vector<esint> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

//	if (!Communication::computeSplitters(_meshData.nIDs, permutation, _nDistribution)) {
//		eslog::error("MESIO internal error: cannot compute permuted nodes splitters.\n");
//	}

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(3 * info::mpi::size + _meshData.nIDs.size() + _meshData.coordinates.size() * sizeof(Point) / sizeof(esint));

	size_t prevsize;
	auto nbegin = permutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of coordinates

		auto n = nbegin;
		for ( ; n != permutation.end() && _meshData.nIDs[*n] < _nDistribution[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n + 1));
		}
		sBuffer[prevsize + 2] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("distribute permuted nodes.\n");
	}

	_meshData.nIDs.clear();
	_meshData.coordinates.clear();

	size_t offset = 0;
	Point point;
	for (int r = 0; r < info::mpi::size; r++) {
		++offset;
		size_t csize = rBuffer[++offset]; // coordinates
		++offset;

		for (size_t c = 0; c < csize; ++c) {
			_meshData.nIDs.push_back(rBuffer[offset++]);
			memcpy(reinterpret_cast<void*>(&point), rBuffer.data() + offset, sizeof(Point));
			_meshData.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(esint);
		}
	}
}

void Input::balanceElements()
{
	std::vector<esint> eCurrent = Communication::getDistribution<esint>(_meshData.esize.size());
	_eDistribution = tarray<esint>::distribute(info::mpi::size, eCurrent.back());

	std::vector<esint> nCurrent = Communication::getDistribution<esint>(_meshData.enodes.size());
	std::vector<esint> nTarget;

	if (info::mpi::rank == 0) {
		nTarget.push_back(0);
	}

	esint nodeOffset = nCurrent[info::mpi::rank];
	size_t eTargetIndex = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), eCurrent[info::mpi::rank] + 1) - _eDistribution.begin();
	for (size_t n = 0; n < _meshData.esize.size(); ++n) {
		nodeOffset += _meshData.esize[n];
		if (eCurrent[info::mpi::rank] + (esint)n + 1 == _eDistribution[eTargetIndex]) {
			nTarget.push_back(nodeOffset);
			++eTargetIndex;
		}
	}
	Communication::allGatherUnknownSize(nTarget);
	nTarget.resize(info::mpi::size + 1, nTarget.back());

	if (!Communication::balance(_meshData.enodes, nCurrent, nTarget)) {
		eslog::internalFailure("balance element nodes.\n");
	}
	if (!Communication::balance(_meshData.esize, eCurrent, _eDistribution)) {
		eslog::internalFailure("balance element sizes.\n");
	}
	if (!Communication::balance(_meshData.eIDs, eCurrent, _eDistribution)) {
		eslog::internalFailure("balance element IDs.\n");
	}
	if (!Communication::balance(_meshData.body, eCurrent, _eDistribution)) {
		eslog::internalFailure("balance element bodies.\n");
	}
	if (!Communication::balance(_meshData.etype, eCurrent, _eDistribution)) {
		eslog::internalFailure("balance element types.\n");
	}
	if (!Communication::balance(_meshData.material, eCurrent, _eDistribution)) {
		eslog::internalFailure("balance element materials.\n");
	}

	auto back = _eDistribution.back();
	if (_meshData.eIDs.size()) {
		back = _meshData.eIDs.back();
	}
	Communication::allGather(&back, _eDistribution.data() + 1, sizeof(back), MPI_BYTE);
	for (size_t i = 1; i < _eDistribution.size(); i++) {
		++_eDistribution[i];
	}
}

void Input::balancePermutedElements()
{
	std::vector<esint> permutation(_meshData.eIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });

//	if (!Communication::computeSplitters(_meshData.eIDs, permutation, _eDistribution)) {
//		eslog::error("MESIO internal error: cannot compute permuted elements splitters.\n");
//	}

	std::vector<esint> edist({ 0 });
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<esint> sBuffer, rBuffer;
	// head, esize, eID, etype, body, material, nodes
	sBuffer.reserve(4 * info::mpi::size + 5 * _meshData.esize.size() + _meshData.enodes.size());

	size_t prevsize;
	auto ebegin = permutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of elements
		sBuffer.push_back(0); // number of elements nodes

		auto e = ebegin;
		for ( ; e != permutation.end() && _meshData.eIDs[*e] < _eDistribution[r + 1]; ++e) {
			sBuffer.push_back(_meshData.esize[*e]);
			sBuffer.push_back(_meshData.eIDs[*e]);
			sBuffer.push_back(_meshData.etype[*e]);
			sBuffer.push_back(_meshData.body[*e]);
			sBuffer.push_back(_meshData.material[*e]);
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + edist[*e], _meshData.enodes.begin() + edist[*e + 1]);
			sBuffer[prevsize + 3] += edist[*e + 1] - edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("distribute permuted elements.\n");
	}

	_meshData.esize.clear();
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.body.clear();
	_meshData.material.clear();
	_meshData.enodes.clear();

	size_t offset = 0;
	for (int r = 0; r < info::mpi::size; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		++offset;
		++offset;

		for (size_t e = 0; e < esize; ++e) {
			_meshData.esize.push_back(rBuffer[offset++]);
			_meshData.eIDs.push_back(rBuffer[offset++]);
			_meshData.etype.push_back(rBuffer[offset++]);
			_meshData.body.push_back(rBuffer[offset++]);
			_meshData.material.push_back(rBuffer[offset++]);
			_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _meshData.esize.back());
			offset += _meshData.esize.back();
		}
	}
}

void Input::sortNodes(bool withElementNodes)
{
	profiler::syncstart("sort_nodes");
	profiler::syncparam("size", _meshData.nIDs.size());
	profiler::syncparam("with_enodes", (int)withElementNodes);
	if (std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end())) {
		profiler::syncend("sort_nodes");
		return;
	}

	std::vector<esint> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	std::sort(_meshData.nIDs.begin(), _meshData.nIDs.end());
	utils::permute(_meshData.coordinates, permutation);
	utils::permute(_nregions, permutation, _nregsize);

	if (_meshData._nranks.size()) {
		std::vector<esint> npermutation(_meshData._nranks.size());
		std::vector<esint> ndist = _meshData._nrankdist;
		for (size_t i = 0, index = 0; i < permutation.size(); i++) {
			_meshData._nrankdist[i + 1] = _meshData._nrankdist[i] + ndist[permutation[i] + 1] - ndist[permutation[i]];
			for (esint n = 0; n < ndist[permutation[i] + 1] - ndist[permutation[i]]; ++n, ++index) {
				npermutation[index] = ndist[permutation[i]] + n;
			}
		}

		utils::permute(_meshData._nranks, npermutation);
	}


	if (withElementNodes) {
		for (size_t n = 0; n < _meshData.enodes.size(); n++) {
			_meshData.enodes[n] = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]) - _meshData.nIDs.begin();
		}
		// not accept IDs that do not start from 0
//		std::vector<esint> backpermutation(_meshData.nIDs.size());
//		std::iota(backpermutation.begin(), backpermutation.end(), 0);
//		std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });
//		for (size_t n = 0; n < _meshData.enodes.size(); n++) {
//			_meshData.enodes[n] = backpermutation[_meshData.enodes[n]];
//		}
	}
	profiler::syncend("sort_nodes");
}

void Input::sortElements()
{
	// this function also removes unsupported elements as their type = 0 = Point
	auto ecomp = [&] (esint i, esint j) {
		if (Mesh::edata[_meshData.etype[i]].type != Mesh::edata[_meshData.etype[j]].type) {
			return static_cast<int>(Mesh::edata[_meshData.etype[i]].type) > static_cast<int>(Mesh::edata[_meshData.etype[j]].type);
		} else {
			return _meshData.eIDs[i] < _meshData.eIDs[j];
		}
	};
	profiler::syncstart("sort_elements");
	profiler::syncparam("size", _meshData.eIDs.size());

	std::vector<esint> permutation(_meshData.eIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	if (!std::is_sorted(permutation.begin(), permutation.end(), ecomp)) {
		std::sort(permutation.begin(), permutation.end(), ecomp);
		sortElements(permutation);
	}

	for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
		_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
			return static_cast<int>(Mesh::edata[e].type) >= type; }) - _meshData.etype.begin()
		);
	}
	profiler::syncend("sort_elements");
}

void Input::sortElements(const std::vector<esint> &permutation)
{
	profiler::syncstart("permute_elements");
	profiler::syncparam("size", permutation.size());
	std::vector<esint> edist;
	if (_meshData._edist.size()) {
		edist.swap(_meshData._edist);
	} else {
		edist = std::vector<esint>({ 0 });
		edist.reserve(_meshData.eIDs.size() + 1);
		for (size_t e = 0; e < _meshData.eIDs.size(); e++) {
			edist.push_back(edist.back() + _meshData.esize[e]);
		}
	}

	utils::permute(_meshData.eIDs, permutation);
	utils::permute(_meshData.esize, permutation);
	utils::permute(_meshData.body, permutation);
	utils::permute(_meshData.etype, permutation);
	utils::permute(_meshData.material, permutation);
	utils::permute(_eregions, permutation, _eregsize);

	std::vector<esint> npermutation(_meshData.enodes.size());
	for (size_t i = 0, index = 0; i < permutation.size(); i++) {
		for (esint n = 0; n < _meshData.esize[i]; ++n, ++index) {
			npermutation[index] = edist[permutation[i]] + n;
		}
	}

	utils::permute(_meshData.enodes, npermutation);
	profiler::syncend("permute_elements");
}


void Input::assignRegions(
		std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs,
		std::vector<esint> &distribution,
		size_t &rsize, std::vector<esint> &rbits)
{
	profiler::syncstart("assign_regions");
	profiler::syncparam("size", IDs.size());
	profiler::syncparam("regions", regions.size());
	rsize = regions.size() / (8 * sizeof(esint)) + 1;
	rbits.resize(rsize * IDs.size());

	std::vector<esint> sBuffer, rBuffer;
	std::vector<std::vector<esint>::const_iterator> rend;
	for (auto region = regions.begin(); region != regions.end(); ++region) {
		rend.push_back(region->second.cbegin());
	}

	for (int t = 0; t < info::mpi::size; t++) {
		size_t prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(t); // target

		size_t r = 0;
		for (auto region = regions.begin(); region != regions.end(); ++region, ++r) {
			size_t reqsize = sBuffer.size();
			sBuffer.push_back(0); // region data size

			auto rprev = rend[r], rbegin = rend[r];
			while (rend[r] < region->second.cend() && *rend[r] < distribution[t + 1]) { ++rend[r]; }
			for (auto id = rbegin; id != rend[r]; rprev = id++) {
				if (*rprev + 1 == *id) { // decode interval as <begin;-end>
					if (sBuffer.back() >= 0) {
						sBuffer.push_back(*id);
					}
					sBuffer.back() = -*id;
				} else {
					sBuffer.push_back(*id);
				}
			}
			sBuffer[reqsize] = sBuffer.size() - reqsize - 1;
		}
		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("assign regions.\n");
	}
	profiler::synccheckpoint("exchange");

	size_t offset = 0;
	for (int t = 0; t < info::mpi::size; t++) {
		++offset; // total size
		++offset; // target

		size_t r = 0;
		for (auto region = regions.begin(); region != regions.end(); ++region, ++r) {
			esint byte = r / (8 * sizeof(esint));
			esint bit = (esint)1 << (r % (8 * sizeof(esint)));

			if (rBuffer[offset]) {
				auto id = std::lower_bound(IDs.begin(), IDs.end(), rBuffer[offset + 1]);
				for (esint i = 0; i < rBuffer[offset]; ++i) {
					if (rBuffer[offset + 1 + i] < 0) { // end of interval
						for (esint ii = rBuffer[offset + 1 + i - 1]; ii <= -rBuffer[offset + 1 + i]; ++ii, ++id) {
							rbits[rsize * (id - IDs.begin()) + byte] |= bit;
						}
					} else {
						while (rBuffer[offset + 1 + i] != *id) { ++id; }
						rbits[rsize * (id - IDs.begin()) + byte] |= bit;
					}
				}
			}
			offset += rBuffer[offset] + 1;
		}
	}

	for (auto region = regions.begin(); region != regions.end(); ++region) {
		region->second.clear();
	}
	profiler::synccheckpoint("rbuffer");
	profiler::syncend("assign_regions");
}

void Input::fillRegions(std::map<std::string, std::vector<esint> > &regions, size_t &rsize, std::vector<esint> &rbits)
{
	profiler::syncstart("fill_regions");
	profiler::syncparam("regions", regions.size());
	profiler::syncparam("size", rbits.size());
	size_t r = 0;
	for (auto region = regions.begin(); region != regions.end(); ++region, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));

		region->second.clear();
		for (size_t i = 0; i < rbits.size() / rsize; ++i) {
			if (rbits[rsize * i + byte] & bit) {
				region->second.push_back(i);
			}
		}
	}
	profiler::syncend("fill_regions");
}

void Input::fillNodes()
{
	profiler::syncstart("fill_nodes");
	profiler::syncparam("size", _meshData.coordinates.size());
	size_t threads = info::env::OMP_NUM_THREADS;

	info::mesh->nodes->size = _meshData.coordinates.size();
	info::mesh->nodes->distribution = tarray<size_t>::distribute(threads, _meshData.coordinates.size());

	info::mesh->nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(info::mesh->nodes->distribution, _meshData.nIDs));
	info::mesh->nodes->coordinates = new serializededata<esint, Point >(1, tarray<Point>(info::mesh->nodes->distribution, _meshData.coordinates.cbegin(), _meshData.coordinates.cend()));

	std::vector<size_t> rdistribution = info::mesh->nodes->distribution, rdatadistribution = info::mesh->nodes->distribution;
	for (size_t t = 1; t < threads; t++) {
		++rdistribution[t];
		if (rdistribution[t] < _meshData._nrankdist.size()) {
			rdatadistribution[t] = _meshData._nrankdist[rdistribution[t]];
		} else {
			rdatadistribution[t] = _meshData._nrankdist[rdistribution[threads] - 1];
		}
	}
	++rdistribution[threads];
	rdatadistribution[threads] = _meshData._nrankdist[rdistribution[threads] - 1];

	info::mesh->nodes->ranks = new serializededata<esint, int>(tarray<esint>(rdistribution, _meshData._nrankdist), tarray<int>(rdatadistribution, _meshData._nranks));

	info::mesh->boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES"));
	info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, _meshData.nIDs.size()));
	std::iota(info::mesh->boundaryRegions.back()->nodes->datatarray().begin(), info::mesh->boundaryRegions.back()->nodes->datatarray().end(), 0);
	profiler::syncend("fill_nodes");
}

void Input::fillElements()
{
	profiler::syncstart("fill_elements");
	size_t estart = info::mesh->dimension == 3 ? 0 : 1;

	if (!_etypeDistribution.size()) {
		for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
			_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
				return static_cast<int>(Mesh::edata[e].type) >= type; }) - _meshData.etype.begin()
			);
		}
	}
	profiler::syncparam("size", _etypeDistribution[estart]);
	_eDistribution = Communication::getDistribution(_etypeDistribution[estart]);

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads), eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, _etypeDistribution[estart]);

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), _meshData._edist.begin() + edistribution[t], _meshData._edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), _meshData._edist.begin() + edistribution[t] + 1, _meshData._edist.begin() + edistribution[t + 1] + 1);
		}

		// till now, IDs are irelevant
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[info::mpi::rank] + edistribution[t]);

		eBody[t].insert(eBody[t].end(), _meshData.body.begin() + edistribution[t], _meshData.body.begin() + edistribution[t + 1]);
		tnodes[t].insert(tnodes[t].end(), _meshData.enodes.begin() + _meshData._edist[edistribution[t]], _meshData.enodes.begin() + _meshData._edist[edistribution[t + 1]]);
		eMat[t].insert(eMat[t].end(), _meshData.material.begin() + edistribution[t], _meshData.material.begin() + edistribution[t + 1]);

		epointers[t].resize(edistribution[t + 1] - edistribution[t]);
		for (size_t e = edistribution[t], i = 0; e < edistribution[t + 1]; ++e, ++i) {
			epointers[t][i] = &Mesh::edata[_meshData.etype[e]];
		}
	}

	info::mesh->elements->distribution.process.offset = _eDistribution[info::mpi::rank];
	info::mesh->elements->distribution.process.next = _eDistribution[info::mpi::rank + 1];
	info::mesh->elements->distribution.process.size = _etypeDistribution[estart];
	info::mesh->elements->distribution.process.totalSize = _etypeDistribution.back();
	info::mesh->elements->distribution.threads = edistribution;
	info::mesh->elements->IDs = new serializededata<esint, esint>(1, eIDs);
	info::mesh->elements->nodes = new serializededata<esint, esint>(tedist, tnodes);
	info::mesh->elements->epointers = new serializededata<esint, Element*>(1, epointers);
	info::mesh->elements->material = new serializededata<esint, int>(1, eMat);
	info::mesh->elements->body = new serializededata<esint, int>(1, eBody);

	info::mesh->elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	info::mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(1, tarray<esint>(threads, info::mesh->elements->distribution.process.size));
	std::iota(info::mesh->elementsRegions.back()->elements->datatarray().begin(), info::mesh->elementsRegions.back()->elements->datatarray().end(), 0);
	profiler::syncend("fill_elements");
}

void Input::fillNeighbors()
{
	std::vector<int> realnranks = _meshData._nranks;
	utils::sortAndRemoveDuplicates(realnranks);

	info::mesh->neighborsWithMe.clear();
	info::mesh->neighborsWithMe.insert(info::mesh->neighborsWithMe.end(), realnranks.begin(), realnranks.end());

	info::mesh->neighbors.clear();
	for (size_t n = 0; n < info::mesh->neighborsWithMe.size(); n++) {
		if (info::mesh->neighborsWithMe[n] != info::mpi::rank) {
			info::mesh->neighbors.push_back(info::mesh->neighborsWithMe[n]);
		}
	}
}

void Input::fillBoundaryRegions()
{
	if (info::ecf->input.omit_face_sets) {
		return;
	}
	profiler::syncstart("fill_boundary_regions");
	size_t threads = info::env::OMP_NUM_THREADS;
	size_t estart = info::mesh->dimension == 3 ? 0 : 1;

	std::vector<std::vector<esint> > tedist(threads), tnodes(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	for (int i = estart; i < 2; i++) {
		std::vector<esint> named;
		for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
			if (eregion->second.size() && _etypeDistribution[estart] <= eregion->second.front()) {
				if (_etypeDistribution[i] <= eregion->second.front() && eregion->second.front() < _etypeDistribution[i + 1]) {
					named.insert(named.end(), eregion->second.begin(), eregion->second.end());
				}
			}
		}
		utils::sortAndRemoveDuplicates(named);
		int hasunnamed = 0, add = 0;
		if (named.size() < (size_t)(_etypeDistribution[i + 1] - _etypeDistribution[i])) {
			hasunnamed = 1;
		}
		Communication::allReduce(&hasunnamed, &add, 1, MPI_INT, MPI_SUM);
		if (add) {
			std::vector<esint> &unnamed = i == 1 ? _meshData.eregions["NAMELESS_EDGE_SET"] : _meshData.eregions["NAMELESS_FACE_SET"];
			auto nit = named.begin();
			for (esint index = _etypeDistribution[i]; index < _etypeDistribution[i + 1] && nit != named.end(); ++index) {
				if (*nit != index) {
					unnamed.push_back(index);
				} else {
					++nit;
				}
			}
			size_t prevsize = unnamed.size();
			esint last = named.size() && named.back() > _etypeDistribution[i] ? named.back() + 1 : _etypeDistribution[i];
			unnamed.resize(_etypeDistribution[i + 1] - _etypeDistribution[i] - named.size());
			std::iota(unnamed.begin() + prevsize, unnamed.end(), last);
		}
	}
	profiler::synccheckpoint("check_unnamed");

	std::vector<int> frominterval, add;
	frominterval.reserve(2 * _meshData.eregions.size());

	for (int i = estart; i < 2; i++) {
		for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
			if (eregion->second.size() && _etypeDistribution[i] <= eregion->second.front() && eregion->second.front() < _etypeDistribution[i + 1]) {
				frominterval.push_back(1);
			} else {
				frominterval.push_back(0);
			}
		}
	}
	add.resize(frominterval.size());
	Communication::allReduce(frominterval.data(), add.data(), frominterval.size(), MPI_INT, MPI_SUM);
	for (int i = estart, j = 0; i < 2; i++) {
		for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
			if (add[j++]) {
				info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(eregion->first));
				info::mesh->boundaryRegions.back()->dimension = info::mesh->boundaryRegions.back()->originalDimension = 2 - i;

				std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, eregion->second.size());
				std::vector<esint> eregiondist(eregion->second.size() + 1);
				for (size_t e = 0; e < eregion->second.size(); e++) {
					eregiondist[e + 1] = eregiondist[e] + _meshData.esize[eregion->second[e]];
				}

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					tedist[t].clear();
					if (t == 0) {
						tedist[t].insert(tedist[t].end(), eregiondist.begin() + edistribution[t], eregiondist.begin() + edistribution[t + 1] + 1);
					} else {
						tedist[t].insert(tedist[t].end(), eregiondist.begin() + edistribution[t] + 1, eregiondist.begin() + edistribution[t + 1] + 1);
					}

					tnodes[t].resize(eregiondist[edistribution[t + 1]] - eregiondist[edistribution[t]]);
					for (size_t e = edistribution[t], index = 0; e < edistribution[t + 1]; ++e) {
						for (esint n = 0; n < _meshData.esize[eregion->second[e]]; ++n, ++index) {
							tnodes[t][index] = _meshData.enodes[_meshData._edist[eregion->second[e]] + n];
						}
					}

					epointers[t].resize(edistribution[t + 1] - edistribution[t]);
					for (size_t e = edistribution[t], i = 0; e < edistribution[t + 1]; ++e, ++i) {
						epointers[t][i] = &Mesh::edata[_meshData.etype[eregion->second[e]]];
					}
				}

				info::mesh->boundaryRegions.back()->distribution.threads = edistribution;
				info::mesh->boundaryRegions.back()->elements = new serializededata<esint, esint>(tedist, tnodes);
				info::mesh->boundaryRegions.back()->epointers = new serializededata<esint, Element*>(1, epointers);
			}
		}
	}
	profiler::synccheckpoint("fill");
	profiler::syncend("fill_boundary_regions");
}

void Input::fillNodeRegions()
{
	profiler::syncstart("fill_node_regions");
	size_t threads = info::env::OMP_NUM_THREADS;

	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion) {
		info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(nregion->first));
		info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, { threads, nregion->second });
	}
	profiler::syncend("fill_node_regions");
}

void Input::fillElementRegions()
{
	profiler::syncstart("fill_element_regions");
	size_t threads = info::env::OMP_NUM_THREADS;
	size_t estart = info::mesh->dimension == 3 ? 0 : 1;

	std::vector<int> fromelements(_meshData.eregions.size()), add(_meshData.eregions.size());
	size_t i = 0;
	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion, ++i) {
		if (eregion->second.size() && eregion->second.front() < _etypeDistribution[estart]) {
			fromelements[i] = 1;
		}
	}
	Communication::allReduce(fromelements.data(), add.data(), fromelements.size(), MPI_INT, MPI_SUM);
	i = 0;
	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion, ++i) {
		auto it = std::lower_bound(eregion->second.begin(), eregion->second.end(), _etypeDistribution.back());
		if (it != eregion->second.end()) {
			eregion->second.resize(it - eregion->second.begin());
		}
		if (add[i]) {
			info::mesh->elementsRegions.push_back(new ElementsRegionStore(eregion->first));
			info::mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(1, { threads, eregion->second });
		}
	}
	profiler::syncend("fill_element_regions");
}

void Input::reindexElementNodes()
{
	profiler::syncstart("reindex_enodes");
	int threads = info::env::OMP_NUM_THREADS;
	int estart = info::mesh->dimension == 3 ? 0 : 1;
	profiler::syncparam("size", _etypeDistribution[estart]);

	if (!_etypeDistribution.size()) {
		for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
			_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
				return static_cast<int>(Mesh::edata[e].type) >= type; }) - _meshData.etype.begin()
			);
		}
	}

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	std::vector<esint> edistribution = tarray<esint>::distribute(threads, _etypeDistribution[estart]);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		for (auto n = _meshData.enodes.begin() + _meshData._edist[edistribution[t]]; n != _meshData.enodes.begin() + _meshData._edist[edistribution[t + 1]]; ++n) {
			*n = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), *n) - _meshData.nIDs.begin();
		}
	}
	profiler::syncend("reindex_enodes");
}

void Input::reindexBoundaryNodes()
{
	profiler::syncstart("reindex_boundary_nodes");
	size_t threads = info::env::OMP_NUM_THREADS;

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->originalDimension) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = info::mesh->boundaryRegions[r]->elements->begin(t)->begin(); n != info::mesh->boundaryRegions[r]->elements->end(t)->begin(); ++n) {
					*n = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), *n) - info::mesh->nodes->IDs->datatarray().begin();
				}
			}
		}
	}
	profiler::syncend("reindex_boundary_nodes");
}

void Input::removeDuplicateElements()
{
	profiler::syncstart("remove_duplicated_elements");
	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	size_t estart = info::mesh->dimension == 3 ? 0 : 1;
	profiler::syncparam("size", _etypeDistribution[estart]);

	std::vector<esint> permutation(_etypeDistribution[estart]);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		edata<esint> e1(_meshData.enodes.data(), _meshData._edist[i], _meshData._edist[i + 1]);
		edata<esint> e2(_meshData.enodes.data(), _meshData._edist[j], _meshData._edist[j + 1]);
		if (e1 == e2) {
			return _meshData.eIDs[i] < _meshData.eIDs[j];
		}
		return e1 < e2;
	});

	profiler::synccheckpoint("sort");

	for (size_t e = 1, unique = 0; e < permutation.size(); ++e) {
		edata<esint> e1(_meshData.enodes.data(), _meshData._edist[permutation[unique]], _meshData._edist[permutation[unique] + 1]);
		edata<esint> e2(_meshData.enodes.data(), _meshData._edist[permutation[e]], _meshData._edist[permutation[e] + 1]);
		if (e1 == e2) {
			_meshData._duplicateElements.push_back(MeshBuilder::Duplicate(_meshData.eIDs[permutation[e]], _meshData.eIDs[permutation[unique]], permutation[e], permutation[unique]));
		} else {
			unique = e;
		}
	}

	profiler::synccheckpoint("compute_duplicates");

	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

	if (info::mesh->neighbors.size()) { // duplicated elements have to have all nodes held by other rank
		for (esint e = 0; e < _etypeDistribution[estart]; ++e) {
			esint node = _meshData.enodes[_meshData._edist[e]];
			std::vector<esint> ranks(_meshData._nranks.begin() + _meshData._nrankdist[node], _meshData._nranks.begin() + _meshData._nrankdist[node + 1]);
			for (esint n = _meshData._edist[e] + 1; n < _meshData._edist[e + 1]; ++n) {
				size_t intersection = 0, current = 0;
				for (int r = _meshData._nrankdist[_meshData.enodes[n]]; r < _meshData._nrankdist[_meshData.enodes[n] + 1]; ++r) {
					while (current < ranks.size() && ranks[current] < _meshData._nranks[r]) { ++current; }
					if (current < ranks.size() && ranks[current] == _meshData._nranks[r]) {
						ranks[intersection++] = ranks[current++];
					}
				}
				ranks.resize(intersection);
			}
			esint roffset = 0;
			for (auto r = ranks.begin(); r != ranks.end(); ++r) {
				if (*r < info::mpi::rank) {
					while (info::mesh->neighbors[roffset] < *r) { ++roffset; }
					sBuffer[roffset].push_back(_meshData.esize[e]);
					for (esint n = _meshData._edist[e]; n < _meshData._edist[e + 1]; ++n) {
							sBuffer[roffset].push_back(_meshData.nIDs[_meshData.enodes[n]]);
					}
					for (size_t n = 0; n < _eregsize; ++n) {
						sBuffer[roffset].push_back(_eregions[_eregsize * e + n]);
					}
					sBuffer[roffset].push_back(e); // save offset for removing
				}
			}
		}
	}
	profiler::synccheckpoint("search_potential");

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("cannot exchange duplicated elements.\n");
	}
	profiler::synccheckpoint("exchange");

	sBuffer.clear();
	sBuffer.resize(info::mesh->neighbors.size());

	for (size_t r = 0; r < rBuffer.size(); ++r) {
		size_t i = 0;
		while (i < rBuffer[r].size()) {
			for (esint n = 0; n < rBuffer[r][i]; ++n) {
				rBuffer[r][i + 1 + n] = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), rBuffer[r][i + 1 + n]) - _meshData.nIDs.begin();
			}
			edata<esint> other(rBuffer[r].data(), i + 1, i + 1 + rBuffer[r][i]);
			auto eit = std::lower_bound(permutation.begin(), permutation.end(), other, [&] (esint i, const edata<esint> &other) {
				edata<esint> e(_meshData.enodes.data(), _meshData._edist[i], _meshData._edist[i + 1]);
				return e < other;
			});
			if (eit != permutation.end()) {
				edata<esint> e(_meshData.enodes.data(), _meshData._edist[*eit], _meshData._edist[*eit + 1]);
				if (e == other) {
					for (size_t n = 0; n < _eregsize; ++n) {
						_eregions[*eit * _eregsize + n] |= rBuffer[r][i + 1 + rBuffer[r][i] + n];
					}
					sBuffer[r].push_back(rBuffer[r][i + 1 + rBuffer[r][i] + _eregsize]);
				}
			}
			i += 1 + rBuffer[r][i] + _eregsize + 1;
		}
	}

	profiler::synccheckpoint("sbuffer");

	rBuffer.clear();
	rBuffer.resize(info::mesh->neighbors.size());
	if (!Communication::receiveLowerUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("cannot exchange duplicated elements.\n");
	}
	profiler::synccheckpoint("exchange");

	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); ++i) {
			_meshData._duplicateElements.push_back(MeshBuilder::Duplicate(_meshData.eIDs[rBuffer[r][i]], -1, rBuffer[r][i], -1));
		}
	}
	profiler::synccheckpoint("rbuffer");

	if (_meshData._duplicateElements.size() == 0) {
		profiler::syncend("remove_duplicated_elements");
		return;
	}

	std::sort(_meshData._duplicateElements.begin(), _meshData._duplicateElements.end(), MeshBuilder::Duplicate());

	profiler::synccheckpoint("sort");

	auto remove = _meshData._duplicateElements.begin();
	esint last = 0, lastn = 0;
	for (size_t e = 0; e < _meshData.eIDs.size(); ++e) {
		while (remove != _meshData._duplicateElements.end() && remove->id < _meshData.eIDs[e]) { ++remove; }
		if (remove == _meshData._duplicateElements.end() || remove->id != _meshData.eIDs[e]) {
			_meshData.eIDs[last] = _meshData.eIDs[e];
			_meshData.etype[last] = _meshData.etype[e];
			_meshData.esize[last] = _meshData.esize[e];
			_meshData.body[last] = _meshData.body[e];
			_meshData.material[last] = _meshData.material[e];
			for (esint n = 0; n < _meshData.esize[last]; ++n) {
				_meshData.enodes[lastn + n] = _meshData.enodes[_meshData._edist[e] + n];
			}
			if (_eregions.size()) {
				for (size_t n = 0; n < _eregsize; ++n) {
					_eregions[_eregsize * last + n] = _eregions[_eregsize * e + n];
				}
			}
			lastn += _meshData.esize[last++];
		} else {
			if (_eregions.size() && remove->target != -1) {
				auto tit = std::lower_bound(_meshData.eIDs.begin(), _meshData.eIDs.begin() + last, remove->target);
				if (tit != _meshData.eIDs.begin() + last && *tit == remove->target) {
					size_t target = tit - _meshData.eIDs.begin();
					for (size_t n = 0; n < _eregsize; ++n) {
						_eregions[_eregsize * target + n] |= _eregions[_eregsize * e + n];
					}
				}
			}
		}
	}
	_meshData.eIDs.resize(last);
	_meshData.etype.resize(last);
	_meshData.esize.resize(last);
	_meshData.body.resize(last);
	_meshData.material.resize(last);
	_meshData.enodes.resize(lastn);
	_meshData._edist.clear();
	if (_eregions.size()) {
		_eregions.resize(_eregsize * last);
	}
	for (auto it = _etypeDistribution.begin(); it != _etypeDistribution.end(); ++it) {
		*it -= _meshData._duplicateElements.size();
	}
	profiler::synccheckpoint("remove");
	profiler::syncend("remove_duplicated_elements");
}

void Input::searchDuplicateNodes()
{
	searchDuplicateNodes(_meshData.coordinates, _meshData.nIDs, [&] (esint id, esint target) {
		_meshData._duplicateNodes.push_back(MeshBuilder::Duplicate{ _meshData.nIDs[id], _meshData.nIDs[target], id, target });
	});
	std::sort(_meshData._duplicateNodes.begin(), _meshData._duplicateNodes.end(), MeshBuilder::Duplicate());
}

void Input::coupleDuplicateNodes()
{
	profiler::syncstart("couple_duplicate_nodes");
	for (auto n = _meshData.enodes.begin(); n != _meshData.enodes.end(); ++n) {
		auto projection = std::lower_bound(_meshData._duplicateNodes.begin(), _meshData._duplicateNodes.end(), *n, [] (MeshBuilder::Duplicate &projection, esint nid) {
			return projection.id < nid;
		});
		if (projection != _meshData._duplicateNodes.end() && projection->id == *n) {
			*n = projection->target;
		}
	}
	profiler::syncend("couple_duplicate_nodes");
}

void Input::searchDuplicateNodes(std::vector<_Point<esfloat> > &coordinates, std::vector<esint> &ids, std::function<void(esint id, esint target)> merge)
{
	profiler::syncstart("search_duplicate_nodes");
	profiler::syncparam("size", coordinates.size());

	KDTree<esfloat> tree(coordinates);

	profiler::synccheckpoint("kdtree_build");
	profiler::param("kdtree_levels", tree.levels);

	double eps = info::ecf->input.duplication_tolerance;
	std::vector<esint> duplicate(tree.permutation.size(), -1);

	for (esint i = std::exp2(tree.levels), first = i; i < std::exp2(tree.levels + 1); ++i) {
		esint begin = tree.begin(i);
		esint end = tree.end(i);
		if (begin == end) {
			continue;
		}

		_Point<esfloat> min;
		tree.boxMin(i, min);

		auto check = [&] (esint p, esint begin, esint end) {
			for (auto pp = tree.permutation.cbegin() + begin; pp != tree.permutation.cbegin() + end; ++pp) {
				if (
						coordinates[tree.permutation[p]].x <= coordinates[*pp].x + eps && coordinates[*pp].x - eps <= coordinates[tree.permutation[p]].x &&
						coordinates[tree.permutation[p]].y <= coordinates[*pp].y + eps && coordinates[*pp].y - eps <= coordinates[tree.permutation[p]].y &&
						coordinates[tree.permutation[p]].z <= coordinates[*pp].z + eps && coordinates[*pp].z - eps <= coordinates[tree.permutation[p]].z) {

					if (duplicate[pp - tree.permutation.cbegin()] >= 0) {
						merge(tree.permutation[p], tree.permutation[duplicate[pp - tree.permutation.cbegin()]]);
						duplicate[p] = duplicate[pp - tree.permutation.cbegin()];
					} else {
						merge(tree.permutation[p], *pp);
						duplicate[p] = pp - tree.permutation.cbegin();
					}
					break;
				}
			}
		};

		std::function<void(size_t, size_t, esint)> traverse = [&] (size_t node, size_t max, esint p) {
			if (coordinates[tree.permutation[p]][tree.splitters[node].d] <= tree.splitters[node].value + eps) {
				if (2 * node < tree.splitters.size()) {
					traverse(2 * node, max, p);
				} else {
					if (2 * node < max) {
						check(p, tree.begin(2 * node), tree.end(2 * node));
					}
				}
			}
			if (tree.splitters[node].value - eps <= coordinates[tree.permutation[p]][tree.splitters[node].d]) {
				if (2 * node < tree.splitters.size()) {
					traverse(2 * node + 1, max, p);
				} else {
					if (2 * node + 1 < max) {
						check(p, tree.begin(2 * node + 1), tree.end(2 * node + 1));
					}
				}
			}
		};

		if (tree.splitters.size() > 1) {
			for (auto p = tree.permutation.cbegin() + begin; p != tree.permutation.cbegin() + end; ++p) {
				if ((coordinates[*p].x <= min.x + eps) || (coordinates[*p].y <= min.y + eps) || (coordinates[*p].z <= min.z + eps)) {
					traverse(1, std::max(first, i), p - tree.permutation.cbegin());
				}
			}
		}

		for (auto left = tree.permutation.cbegin() + begin, right = left + 1; right != tree.permutation.cbegin() + end; ++right) {
			if (duplicate[right - tree.permutation.cbegin()] != -1) {
				continue;
			}
			while (left != right && (coordinates[*left].x + eps < coordinates[*right].x || duplicate[left - tree.permutation.cbegin()] != -1)) {
				++left;
			}
			for (auto mid = left; mid != right;) {
				while (mid != right && coordinates[*mid].y + eps <  coordinates[*right].y) {
					++mid;
				}
				if (mid != right && coordinates[*mid].y - eps <= coordinates[*right].y) {
					if (duplicate[mid - tree.permutation.cbegin()] == -1 && coordinates[*right].z <= coordinates[*mid].z + eps && coordinates[*mid].z - eps <= coordinates[*right].z) {
						if (duplicate[mid - tree.permutation.cbegin()] >= 0) {
							merge(*right, tree.permutation[duplicate[mid - tree.permutation.cbegin()]]);
							duplicate[right - tree.permutation.cbegin()] = duplicate[mid - tree.permutation.cbegin()];
						} else {
							merge(*right, *mid);
							duplicate[right - tree.permutation.cbegin()] = mid - tree.permutation.cbegin();
						}
					}
					++mid;
				}
				while (mid != right && coordinates[*right].y + eps < coordinates[*mid].y) {
					++mid;
				}
			}
		}
	}

	profiler::synccheckpoint("merge");
	profiler::syncend("search_duplicate_nodes");
}
