
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

namespace espreso {
namespace builder {

void assignBuckets(const OrderedMesh &mesh, const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered)
{
	clustered.buckets = sfc.buckets(sfc.depth);
	std::vector<size_t> cdistribution = tarray<size_t>::distribute(info::env::OMP_NUM_THREADS, mesh.coordinates.size());
	std::vector<size_t> edistribution = tarray<size_t>::distribute(info::env::OMP_NUM_THREADS, mesh.etype.size());

	clustered.nbuckets.resize(mesh.coordinates.size());
	clustered.ebuckets.resize(mesh.etype.size());

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {

		// nbuckets -> just ask SFC to get bucket
		for (size_t n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			clustered.nbuckets[n] = sfc.getBucket(mesh.coordinates[n]);
		}

		// ebuckets -> we need to choose a node and ask a neighboring process to a bucket
		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			clustered.ebuckets[e] = mesh.enodes[mesh.edist[e]];
			for (esint n = mesh.edist[e]; n < mesh.edist[e + 1]; n++) {
				if (mesh.enodes[n] / mesh.nchunk == info::mpi::rank) {
					clustered.ebuckets[e] = clustered.nbuckets[mesh.enodes[n] - mesh.noffset];
				} else {
					if (std::abs(clustered.ebuckets[e] / mesh.nchunk - info::mpi::rank) > std::abs(mesh.enodes[n] / mesh.nchunk - info::mpi::rank)) {
						clustered.ebuckets[e] = mesh.enodes[n];
					}
				}
			}
		}
	}

	// TODO: measure if it is better to avoid sorting
	std::vector<esint, initless_allocator<esint> > permutation(clustered.ebuckets.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (const esint &i, const esint &j) { return clustered.ebuckets[i] < clustered.ebuckets[j]; });

	size_t skipped = 0;
	std::vector<esint> sBuffer, rBuffer;
	{ // build buffer with required nodes
		sBuffer.reserve(permutation.size() + 3 * info::mpi::size);
		auto p = permutation.begin();
		for (int r = 0; r < info::mpi::size; r++) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0);
			sBuffer.push_back(r);
			sBuffer.push_back(info::mpi::rank);

			if (r == info::mpi::rank) {
				while (p != permutation.end() && clustered.ebuckets[*p] < mesh.nchunk * (r + 1)) {
					++skipped;
					++p;
				}
			} else {
				while (p != permutation.end() && clustered.ebuckets[*p] < mesh.nchunk * (r + 1)) {
					sBuffer.push_back(clustered.ebuckets[*p++]);
				}
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("ask neighbors for nodes buckets.\n");
	}

	{ // requests are not sorted according to rank -> reorder them
		std::vector<esint>().swap(sBuffer);
		sBuffer.resize(rBuffer.size());

		std::vector<esint, initless_allocator<esint> > toRankOffset(info::mpi::size + 1);
		auto buffer = rBuffer.begin();
		for (int r = 0; r < info::mpi::size; ++r) {
			toRankOffset[*(buffer + 2)] = *buffer;
			buffer += *buffer;
		}
		utils::sizesToOffsets(toRankOffset.data(), toRankOffset.data() + toRankOffset.size(), 0);

		for (size_t offset = 0; offset < rBuffer.size(); ) {
			esint size = rBuffer[offset++];
			++offset; // my rank
			int toRank = rBuffer[offset++];
			sBuffer[toRankOffset[toRank] + 0] = size;
			sBuffer[toRankOffset[toRank] + 1] = toRank;
			sBuffer[toRankOffset[toRank] + 2] = info::mpi::rank;
			for (esint n = toRankOffset[toRank] + 3; n < toRankOffset[toRank + 1]; ++n) {
				sBuffer[n] = clustered.nbuckets[rBuffer[offset++] - mesh.noffset];
			}
		}
	}

	// check if a point-to-point version can be used for a better performance
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("return nodes buckets.\n");
	}
	std::vector<esint>().swap(sBuffer);

	{ // returned nodes are in different order than send requests -> reorder them
		std::vector<esint, initless_allocator<esint> > toRankOffset(info::mpi::size + 1);
		auto buffer = rBuffer.begin();
		for (int r = 0; r < info::mpi::size; ++r) {
			toRankOffset[*(buffer + 2)] = *buffer;
			buffer += *buffer;
		}
		utils::sizesToOffsets(toRankOffset.data(), toRankOffset.data() + toRankOffset.size(), 0);

		size_t p = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			esint offset = toRankOffset[r];
			esint size = rBuffer[offset++];
			offset++; // my rank
			offset++; // r

			if (r == info::mpi::rank) {
				p += skipped;
			} else {
				for (esint n = 3; n < size; ++n) {
					clustered.ebuckets[permutation[p++]] = rBuffer[offset++];
				}
			}
		}
	}
}

void clusterize(OrderedMesh &mesh, ClusteredMesh &clustered)
{
	clustered.dimension = mesh.dimension;
	std::vector<esint, initless_allocator<esint> > npermutation(clustered.nbuckets.size()), epermutation(clustered.ebuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (const esint &i, const esint &j) { return clustered.nbuckets[i] < clustered.nbuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (const esint &i, const esint &j) { return clustered.ebuckets[i] < clustered.ebuckets[j]; });

	if (!Communication::computeSplitters(clustered.ebuckets, epermutation, clustered.splitters, mesh.etotal, clustered.buckets)) {
		eslog::internalFailure("cannot balance according to SFC.\n");
	}

	std::vector<esint, initless_allocator<esint> > nborders(info::mpi::size + 1), eborders(info::mpi::size + 1);
	nborders[0] = eborders[0] = 0;
	for (int r = 0; r < info::mpi::size; ++r) {
		nborders[r + 1] = nborders[r];
		while (nborders[r + 1] < mesh.nsize && clustered.nbuckets[npermutation[nborders[r + 1]]] < clustered.splitters[r + 1]) {
			++nborders[r + 1];
		}
		eborders[r + 1] = eborders[r];
		while (eborders[r + 1] < mesh.esize && clustered.ebuckets[epermutation[eborders[r + 1]]] < clustered.splitters[r + 1]) {
			++eborders[r + 1];
		}
	}

	#pragma omp parallel for
	for (int r = 0; r < info::mpi::size; ++r) {
		std::sort(npermutation.begin() + nborders[r], npermutation.begin() + nborders[r + 1]);
		std::sort(epermutation.begin() + eborders[r], epermutation.begin() + eborders[r + 1]);
	}

	std::vector<esint, initless_allocator<esint> > sBuffer, rBuffer;

	{ // compute size of the send buffer
		size_t ssize = 0;
		ssize += 6 * info::mpi::size;
		ssize += mesh.nsize;
		ssize += mesh.esize;
		ssize += mesh.enodes.size();
		for (int r = 0; r < info::mpi::size; ++r) {
			ssize += utils::reinterpret_size<esint, _Point<esfloat> >(nborders[r + 1] - nborders[r]);
			ssize += utils::reinterpret_size<esint, char>(eborders[r + 1] - eborders[r]);
		}
		sBuffer.reserve(ssize);
	}

	for (int r = 0; r < info::mpi::size; ++r) {
		size_t prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(info::mpi::rank);
		sBuffer.push_back(nborders[r + 1] - nborders[r]);
		sBuffer.push_back(eborders[r + 1] - eborders[r]); // number of elements
		sBuffer.push_back(0); // number of elements nodes

		for (esint n = nborders[r]; n < nborders[r + 1]; ++n) {
			sBuffer.push_back(npermutation[n] + mesh.noffset);
		}
		char *pbuffer = reinterpret_cast<char*>(sBuffer.data() + sBuffer.size());
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, _Point<esfloat> >(nborders[r + 1] - nborders[r]));
		for (esint n = nborders[r]; n < nborders[r + 1]; ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, mesh.coordinates.data() + npermutation[n], sizeof(_Point<esfloat>));
		}

		char *tbuffer = reinterpret_cast<char*>(sBuffer.data() + sBuffer.size());
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, char>(eborders[r + 1] - eborders[r]));
		for (esint e = eborders[r]; e < eborders[r + 1]; ++e, ++tbuffer) {
			memcpy(tbuffer, mesh.etype.data() + epermutation[e], sizeof(Element::CODE));
		}
		for (esint e = eborders[r]; e < eborders[r + 1]; ++e) {
			sBuffer.push_back(epermutation[e] + mesh.eoffset);
			for (esint en = mesh.edist[epermutation[e]]; en < mesh.edist[epermutation[e] + 1]; ++en) {
				sBuffer.push_back(mesh.enodes[en]);
			}
			sBuffer[prevsize + 5] += mesh.edist[epermutation[e] + 1] - mesh.edist[epermutation[e]];
		}
		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}
	mesh.clearNodes();
	mesh.clearElements();
	mesh.clearRegions();
	mesh.clearValues();
	utils::clearVectors(mesh.edist, nborders, eborders, npermutation, epermutation);

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("distribute elements according to SFC.\n");
	}
	utils::clearVector(sBuffer);

	std::vector<size_t, initless_allocator<size_t> > roffset(info::mpi::size);
	{ // preallocate vectors
		esint nodesTotal = 0, elementsTotal = 0, enodesTotal = 0;
		for (size_t offset = 0; offset < rBuffer.size();) {
			++offset; // totalsize
			++offset; // my rank
			esint source = rBuffer[offset++];
			roffset[source] = offset;
			esint nodes = rBuffer[offset++];
			esint elements = rBuffer[offset++];
			esint enodes = rBuffer[offset++];
			offset += nodes + utils::reinterpret_size<esint, _Point<esfloat> >(nodes);
			offset += elements + enodes + utils::reinterpret_size<esint, Element::CODE>(elements);
			nodesTotal += nodes;
			elementsTotal += elements;
			enodesTotal += enodes;
		}
		clustered.noffsets.reserve(nodesTotal);
		clustered.coordinates.reserve(nodesTotal);
		clustered.eoffsets.reserve(elementsTotal);
		clustered.etype.reserve(elementsTotal);
		clustered.enodes.reserve(enodesTotal);
		clustered.edist.reserve(elementsTotal + 1);
	}

	clustered.edist.push_back(0);
	for (int r = 0; r < info::mpi::size; ++r) {
		size_t offset = roffset[r];
		esint nodes = rBuffer[offset++];
		esint elements = rBuffer[offset++];
		++offset; // enodes
		clustered.noffsets.insert(clustered.noffsets.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + nodes);
		offset += nodes;
		clustered.coordinates.insert(clustered.coordinates.end(), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset) + nodes);
		offset += utils::reinterpret_size<esint, _Point<esfloat> >(nodes);

		clustered.etype.insert(clustered.etype.end(), reinterpret_cast<Element::CODE*>(rBuffer.data() + offset), reinterpret_cast<Element::CODE*>(rBuffer.data() + offset) + elements);
		offset += utils::reinterpret_size<esint, Element::CODE>(elements);
		for (esint e = 0; e < elements; ++e) {
			clustered.eoffsets.push_back(rBuffer[offset++]);
			for (int en = 0; en < Mesh::edata[(int)clustered.etype[clustered.eoffsets.size() - 1]].nodes; ++en) {
				clustered.enodes.push_back(rBuffer[offset++]);
			}
			clustered.edist.push_back(clustered.enodes.size());
		}
	}
}

void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered)
{
	std::vector<std::pair<size_t, size_t> > neighbors;

	size_t index = clustered.splitters[info::mpi::rank];
	size_t last = clustered.splitters[info::mpi::rank + 1];
	while (index < last) {
		size_t depth = sfc.depth, bsize = 1;
		while (depth > 1 && index % (bsize * sfc.bucketSize()) == 0 && index + (bsize * sfc.bucketSize()) < last) {
			--depth;
			bsize *= sfc.bucketSize();
		}
		sfc.addSFCNeighbors(depth, index / bsize, clustered.splitters, neighbors);
		index += bsize;
	}

	utils::sortAndRemoveDuplicates(neighbors);

	for (size_t i = 0; i < neighbors.size(); i++) {
		size_t bstep = sfc.buckets(sfc.depth) / sfc.buckets(neighbors[i].first);
		neighbors[i].first = neighbors[i].second * bstep;
		neighbors[i].second = neighbors[i].second * bstep + bstep;
	}

	std::sort(neighbors.begin(), neighbors.end());

	// we assume that where is not any interval across processes (assured by the above addNeighbors alg.)
	auto rank = clustered.splitters.begin();
	for (size_t i = 0; i < neighbors.size(); i++) {
		while (*rank <= (esint)neighbors[i].first) { ++rank; }
		int r = rank - clustered.splitters.begin() - 1;
		if (r != info::mpi::rank && (clustered.neighbors.empty() || clustered.neighbors.back() != r)) {
			clustered.neighbors.push_back(r);
		}
	}
}

}
}
