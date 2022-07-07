
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

namespace espreso {
namespace builder {

void assignBuckets(OrderedNodesBalanced &nodes, OrderedElementsBalanced &elements, const HilbertCurve<esfloat> &sfc, ivector<esint> &nbuckets, ivector<esint> &ebuckets)
{
	nbuckets.resize(nodes.size);
	ebuckets.resize(elements.size);

	esint local = 0;
	ivector<esint> &closest = ebuckets; // reuse memory

	// nbuckets -> just ask SFC to get bucket
	for (esint n = 0; n < nodes.size; ++n) {
		nbuckets[n] = sfc.getBucket(nodes.coordinates[n]);
	}

	// ebuckets -> we need to choose a node and ask a neighboring process to a bucket
	for (esint e = 0, eoffset = 0; e < elements.size; eoffset += Mesh::element(elements.etype[e++]).nodes) {
		closest[e] = elements.enodes[eoffset];
		for (esint n = eoffset; n < eoffset + Mesh::element(elements.etype[e]).nodes; n++) {
			if (elements.enodes[n] / nodes.chunk == info::mpi::rank) {
				closest[e] = elements.enodes[n];
				++local;
				break;
			} else {
				if (std::abs(ebuckets[e] / nodes.chunk - info::mpi::rank) > std::abs(elements.enodes[n] / nodes.chunk - info::mpi::rank)) {
					closest[e] = elements.enodes[n];
				}
			}
		}
	}

	// TODO: measure if it is better to avoid sorting
	std::vector<esint, initless_allocator<esint> > permutation(closest.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (const esint &i, const esint &j) { return closest[i] < closest[j]; });

	ivector<esint> sBuffer, rBuffer;
	{ // build buffer with required nodes
		sBuffer.reserve(permutation.size() + 3 * info::mpi::size);
		auto p = permutation.begin();
		for (int r = 0; r < info::mpi::size; r++) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0);
			sBuffer.push_back(r);
			sBuffer.push_back(info::mpi::rank);

			if (r == info::mpi::rank) {
				p += local;
			} else {
				while (p != permutation.end() && closest[*p] < nodes.chunk * (r + 1)) {
					sBuffer.push_back(closest[*p++]);
				}
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("ask neighbors for nodes buckets.\n");
	}

	{ // requests are not sorted according to rank -> reorder them
		utils::clearVector(sBuffer);
		sBuffer.reserve(rBuffer.size());

		ivector<esint> toRankOffset(info::mpi::size);
		auto buffer = rBuffer.begin();
		for (int r = 0; r < info::mpi::size; ++r) {
			toRankOffset[*(buffer + 2)] = buffer - rBuffer.begin();
			buffer += *buffer;
		}

		for (int r = 0; r < info::mpi::size; ++r) {
			sBuffer.push_back(rBuffer[toRankOffset[r]]);
			sBuffer.push_back(r);
			sBuffer.push_back(info::mpi::rank);
			for (esint n = 3; n < rBuffer[toRankOffset[r]]; ++n) {
				sBuffer.push_back(nbuckets[rBuffer[toRankOffset[r] + n] - nodes.offset]);
			}
		}
	}

	// check if a point-to-point version can be used for a better performance
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("return nodes buckets.\n");
	}
	utils::clearVector(sBuffer);

	{ // returned nodes are in different order than send requests -> reorder them
		ivector<esint> toRankOffset(info::mpi::size);
		auto buffer = rBuffer.begin();
		for (int r = 0; r < info::mpi::size; ++r) {
			toRankOffset[*(buffer + 2)] = buffer - rBuffer.begin();
			buffer += *buffer;
		}

		size_t p = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			esint offset = toRankOffset[r];
			esint size = rBuffer[offset++];
			++offset; // my rank
			++offset; // r

			if (r == info::mpi::rank) {
				for (esint i = 0; i < local; ++i, ++p) {
					ebuckets[permutation[p]] = nbuckets[closest[permutation[p]] - nodes.offset];
				}
			} else {
				for (esint n = 3; n < size; ++n) {
					ebuckets[permutation[p++]] = rBuffer[offset++];
				}
			}
		}
	}
}

void clusterize(OrderedNodesBalanced &inNodes, OrderedElementsBalanced &inElements, ivector<esint> &nbuckets, ivector<esint> &ebuckets, esint buckets, ClusteredNodes &outNodes, ClusteredElements &outElements, ivector<esint> &splitters)
{
	std::vector<esint, initless_allocator<esint> > npermutation(nbuckets.size()), epermutation(ebuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (const esint &i, const esint &j) { return nbuckets[i] < nbuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (const esint &i, const esint &j) { return ebuckets[i] < ebuckets[j]; });

	if (!Communication::computeSplitters(ebuckets, epermutation, splitters, inElements.total, buckets)) {
		eslog::internalFailure("cannot balance according to SFC.\n");
	}

	std::vector<esint, initless_allocator<esint> > nborders(info::mpi::size + 1), eborders(info::mpi::size + 1);
	nborders[0] = eborders[0] = 0;
	for (int r = 0; r < info::mpi::size; ++r) {
		nborders[r + 1] = nborders[r];
		while (nborders[r + 1] < inNodes.size && nbuckets[npermutation[nborders[r + 1]]] < splitters[r + 1]) {
			++nborders[r + 1];
		}
		eborders[r + 1] = eborders[r];
		while (eborders[r + 1] < inElements.size && ebuckets[epermutation[eborders[r + 1]]] < splitters[r + 1]) {
			++eborders[r + 1];
		}
	}

	#pragma omp parallel for
	for (int r = 0; r < info::mpi::size; ++r) {
		std::sort(npermutation.begin() + nborders[r], npermutation.begin() + nborders[r + 1]);
		std::sort(epermutation.begin() + eborders[r], epermutation.begin() + eborders[r + 1]);
	}

	ivector<esint> edist(inElements.size + 1);
	edist[0] = 0;
	for (esint e = 0; e < inElements.size; ++e) {
		edist[e + 1] = edist[e] + Mesh::element(inElements.etype[e]).nodes;
	}

	ivector<esint> sBuffer, rBuffer;

	{ // compute size of the send buffer
		size_t ssize = 0;
		ssize += 6 * info::mpi::size;
		ssize += inNodes.size;
		ssize += inElements.size;
		ssize += inElements.enodes.size();
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
			sBuffer.push_back(npermutation[n] + inNodes.offset);
		}
		char *pbuffer = reinterpret_cast<char*>(sBuffer.data() + sBuffer.size());
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, _Point<esfloat> >(nborders[r + 1] - nborders[r]));
		for (esint n = nborders[r]; n < nborders[r + 1]; ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, inNodes.coordinates.data() + npermutation[n], sizeof(_Point<esfloat>));
		}

		char *tbuffer = reinterpret_cast<char*>(sBuffer.data() + sBuffer.size());
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, char>(eborders[r + 1] - eborders[r]));
		for (esint e = eborders[r]; e < eborders[r + 1]; ++e, ++tbuffer) {
			memcpy(tbuffer, inElements.etype.data() + epermutation[e], sizeof(Element::CODE));
		}
		for (esint e = eborders[r]; e < eborders[r + 1]; ++e) {
			sBuffer.push_back(epermutation[e] + inElements.offset);
			for (esint en = edist[epermutation[e]]; en < edist[epermutation[e] + 1]; ++en) {
				sBuffer.push_back(inElements.enodes[en]);
			}
			sBuffer[prevsize + 5] += edist[epermutation[e] + 1] - edist[epermutation[e]];
		}
		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	utils::clearVectors(inNodes.coordinates);
	utils::clearVectors(inElements.etype, inElements.enodes);
	utils::clearVectors(nborders, eborders, npermutation, epermutation);

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("distribute elements according to SFC.\n");
	}
	utils::clearVector(sBuffer);

	ivector<size_t> roffset(info::mpi::size);
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
		// another space for shared boundary elements that can be added later
		elementsTotal += 100;
		enodesTotal += 8 * 100;

		outNodes.offsets.reserve(nodesTotal);
		outNodes.coordinates.reserve(nodesTotal);
		outElements.offsets.reserve(elementsTotal);
		outElements.etype.reserve(elementsTotal);
		outElements.enodes.reserve(enodesTotal);
	}

	for (int r = 0; r < info::mpi::size; ++r) {
		size_t offset = roffset[r];
		esint nodes = rBuffer[offset++];
		esint elements = rBuffer[offset++];
		++offset; // enodes
		outNodes.offsets.insert(outNodes.offsets.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + nodes);
		offset += nodes;
		outNodes.coordinates.insert(outNodes.coordinates.end(), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset) + nodes);
		offset += utils::reinterpret_size<esint, _Point<esfloat> >(nodes);

		outElements.etype.insert(outElements.etype.end(), reinterpret_cast<Element::CODE*>(rBuffer.data() + offset), reinterpret_cast<Element::CODE*>(rBuffer.data() + offset) + elements);
		offset += utils::reinterpret_size<esint, Element::CODE>(elements);
		for (esint e = 0; e < elements; ++e) {
			outElements.offsets.push_back(rBuffer[offset++]);
			for (int en = 0; en < Mesh::element(outElements.etype[outElements.offsets.size() - 1]).nodes; ++en) {
				outElements.enodes.push_back(rBuffer[offset++]);
			}
		}
	}
}

void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, std::vector<int> &sfcNeighbors)
{
	std::vector<std::pair<size_t, size_t> > neighbors;

	size_t index = splitters[info::mpi::rank];
	size_t last = splitters[info::mpi::rank + 1];
	while (index < last) {
		size_t depth = sfc.depth, bsize = 1;
		while (depth > 1 && index % (bsize * sfc.bucketSize()) == 0 && index + (bsize * sfc.bucketSize()) < last) {
			--depth;
			bsize *= sfc.bucketSize();
		}
		sfc.addSFCNeighbors(depth, index / bsize, splitters, neighbors);
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
	auto rank = splitters.begin();
	for (size_t i = 0; i < neighbors.size(); i++) {
		while (*rank <= (esint)neighbors[i].first) { ++rank; }
		int r = rank - splitters.begin() - 1;
		if (r != info::mpi::rank && (sfcNeighbors.empty() || sfcNeighbors.back() != r)) {
			sfcNeighbors.push_back(r);
		}
	}

	size_t maxNeigh = sfcNeighbors.size();
	Communication::allReduce(&maxNeigh, nullptr, 1, MPITools::getType(maxNeigh).mpitype, MPI_MAX);
	eslog::info(" == MAX SFC NEIGHBORS %69lu == \n", maxNeigh);
}

}
}
