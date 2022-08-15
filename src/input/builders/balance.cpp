
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

static bool chunk(const esint &mpichunk, const int &rank, const std::vector<DatabaseOffset> &blocks, std::vector<DatabaseOffset>::const_iterator &it, esint &begin, esint &end)
{
	if (it != blocks.end() && end == it->local + it->size) {
		++it;
	}
	if (it == blocks.end()) {
		return false;
	}
	if (it->global + end - it->local == mpichunk * (rank + 1)) {
		return false;
	}
	begin = end = 0;
	while (it != blocks.end() && it->global + it->size < mpichunk * rank) { ++it; }
	if (it != blocks.end() && it->global < mpichunk * (rank + 1)) {
		begin = std::max(it->global, mpichunk * rank);
		end = std::min(it->global + it->size, mpichunk * (rank + 1));
		begin = it->local + begin - it->global;
		end = it->local + end - it->global;
	}
	return begin != end;
}

static void distribute(BalancedDistribution &distribution, const esint &total)
{
	distribution.total = total;
	distribution.size = distribution.chunk = total / info::mpi::size + ((total % info::mpi::size) ? 1 : 0);
	distribution.offset = std::min(distribution.chunk * info::mpi::rank, distribution.total);
	if (distribution.total <= distribution.offset + distribution.size) {
		distribution.size = distribution.total - distribution.offset;
	}
}

void balanceFEM(NodesBlocks &inNodes, ElementsBlocks &inElements, OrderedNodesBalanced &outNodes, OrderedElementsBalanced &outElements)
{
	esint total[2] = { (esint)inNodes.coordinates.size(), (esint)inElements.etype.size() };
	Communication::allReduce(total, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);

	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", total[0]);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", total[1]);

	distribute(outNodes, total[0]);
	distribute(outElements, total[1]);

	std::vector<esint> sBuffer, rBuffer, edist;
	edist.reserve(inElements.etype.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < inElements.etype.size(); ++e) {
		edist.push_back(edist.back() + Element::encode(inElements.etype[e]).nodes);
	}

	auto blockcomp = [] (const DatabaseOffset &me, const DatabaseOffset &other) { return me.global < other.global; };
	std::sort(inNodes.blocks.begin(), inNodes.blocks.end(), blockcomp);
	std::sort(inElements.blocks.begin(), inElements.blocks.end(), blockcomp);

	{ // compute size of the send buffer
		size_t ssize = 0;
		std::vector<DatabaseOffset>::const_iterator nit = inNodes.blocks.cbegin();
		std::vector<DatabaseOffset>::const_iterator eit = inElements.blocks.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			ssize += 4;
			while (chunk(outNodes.chunk, r, inNodes.blocks, nit, nbegin, nend)) {
				ssize += 2;
				ssize += utils::reinterpret_size<esint, _Point<esfloat> >(nend - nbegin);
			}
			while (chunk(outElements.chunk, r, inElements.blocks, eit, ebegin, eend)) {
				ssize += 3;
				ssize += utils::reinterpret_size<esint, Element::CODE>(eend - ebegin);
				ssize += edist[eend] - edist[ebegin];
			}
		}
		sBuffer.reserve(ssize);
	}

	{ // build the send buffer
		std::vector<DatabaseOffset>::const_iterator nit = inNodes.blocks.cbegin();
		std::vector<DatabaseOffset>::const_iterator eit = inElements.blocks.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0); // total size
			sBuffer.push_back(r); // target
			sBuffer.push_back(0); // nodes
			sBuffer.push_back(0); // elements

			while (chunk(outNodes.chunk, r, inNodes.blocks, nit, nbegin, nend)) {
				sBuffer.push_back(nit->global + nbegin - nit->local);
				sBuffer.push_back(nend - nbegin);
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(inNodes.coordinates.data() + nbegin), utils::reinterpret_end<esint>(inNodes.coordinates.data() + nbegin, nend - nbegin));
				++sBuffer[prevsize + 2];
			}
			while (chunk(outElements.chunk, r, inElements.blocks, eit, ebegin, eend)) {
				sBuffer.push_back(eit->global + ebegin - eit->local);
				sBuffer.push_back(eend - ebegin);
				sBuffer.push_back(edist[eend] - edist[ebegin]);
				// it can cause Invalid read message with Valgrind (we touch after etype array but we never use that data)
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(inElements.etype.data() + ebegin), utils::reinterpret_end<esint>(inElements.etype.data() + ebegin, eend - ebegin));
				sBuffer.insert(sBuffer.end(), inElements.enodes.data() + edist[ebegin], inElements.enodes.data() + edist[eend]);
				++sBuffer[prevsize + 3];
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}

		utils::clearVector(edist);
		utils::clearVectors(inNodes.coordinates, inNodes.blocks);
		utils::clearVectors(inElements.etype, inElements.enodes, inElements.blocks);
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot balance parsed data.\n");
	}
	utils::clearVector(sBuffer);

	outNodes.coordinates.resize(outNodes.size);
	outElements.etype.resize(outElements.size);

	size_t offset = 0, enodesTotal = 0;
	std::vector<std::pair<esint, esint> > eintervals; // element offset, rBuffer offset
	for (int r = 0; r < info::mpi::size; ++r) {
		++offset; // totalsize
		++offset; // my rank
		esint nblocks = rBuffer[offset++];
		esint eblocks = rBuffer[offset++];
		for (esint n = 0; n < nblocks; ++n) {
			esint coffset = rBuffer[offset++];
			esint csize = rBuffer[offset++];
			memcpy(outNodes.coordinates.data() + coffset - outNodes.offset, reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset), csize * sizeof(_Point<esfloat>));
			offset += utils::reinterpret_size<esint, _Point<esfloat> >(csize);
		}
		for (esint e = 0; e < eblocks; ++e) {
			esint eoffset = rBuffer[offset++];
			eintervals.push_back(std::make_pair(eoffset, offset));
			esint esize = rBuffer[offset++];
			esint enodes = rBuffer[offset++];
			offset += utils::reinterpret_size<esint, Element::CODE>(esize) + enodes;
			enodesTotal += enodes;
		}
	}
	outElements.enodes.reserve(enodesTotal);
	std::sort(eintervals.begin(), eintervals.end());
	for (size_t i = 0; i < eintervals.size(); ++i) {
		esint eoffset = eintervals[i].first;
		esint esize = rBuffer[eintervals[i].second++];
		eintervals[i].second++; // enodes
		memcpy(outElements.etype.data() + eoffset - outElements.offset, reinterpret_cast<Element::CODE*>(rBuffer.data() + eintervals[i].second), esize * sizeof(Element::CODE));
		eintervals[i].second += utils::reinterpret_size<esint, Element::CODE>(esize);
		for (esint en = 0; en < esize; ++en) {
			for (int nn = 0; nn < Element::encode(outElements.etype[eoffset - outElements.offset + en]).nodes; ++nn) {
				outElements.enodes.push_back(rBuffer[eintervals[i].second++]);
			}
		}
	}
}

template <typename T>
esint insert(esint offset, esint size, esint rbegin, esint rend, T* data, ivector<esint> &buffer)
{
	// we start from beginning, hence, 'offset' is never before 'rbegin'
	if (rend <= offset) {
		return 0;
	}
	T *begin = data, *end = begin + std::min(size, rend - offset);
	buffer.insert(buffer.end(), reinterpret_cast<esint*>(begin), utils::reinterpret_end<esint>(data, end - begin));
	return end - begin;
}

void balanceFVM(NodesBlocks &inNodes, FacesBlocks &inFaces, OrderedNodesBalanced &outNodes, OrderedFacesBalanced &outFaces)
{
	std::vector<esint> sum(4), offset = { (esint)inNodes.coordinates.size(), (esint)inFaces.ftype.size(), (esint)inFaces.owner.size(), (esint)inFaces.neighbor.size() };
	Communication::exscan(sum, offset);

	std::vector<esint> fdistribution = Communication::getDistribution<esint>(sum[1]);
	distribute(outNodes, sum[0]);
	distribute(outFaces.elements, inFaces.elements.blocks.back().size);

	// 1. balance description of faces
	ivector<esint> sBuffer, rBuffer, edist;
	edist.reserve(inFaces.ftype.size() + 2);
	edist.push_back(0);
	for (size_t e = 0; e < inFaces.ftype.size(); ++e) {
		edist.push_back(edist.back() + Element::encode(inFaces.ftype[e]).nodes);
	}
	edist.push_back(edist.back());

	sBuffer.reserve(
			10 * info::mpi::size + // 8 + 2 * reinterpter_size rounding
			utils::reinterpret_size<esint, _Point<esfloat> >(inNodes.coordinates.size()) + utils::reinterpret_size<esint, Element::CODE>(inFaces.ftype.size()) +
			inFaces.fnodes.size() + inFaces.owner.size() + inFaces.neighbor.size());
	for (esint r = 0, coffset = 0, foffset = 0, ooffset = 0, noffset = 0; r < info::mpi::size; ++r) {
		size_t size = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(info::mpi::rank); // source
		sBuffer.push_back(0); // coordinates
		sBuffer.push_back(0); // faces
		sBuffer.push_back(0); // nodes
		sBuffer.push_back(0); // owners
		sBuffer.push_back(0); // neighbors

		sBuffer[size + 3] = insert(offset[0] + coffset, inNodes.coordinates.size() - coffset, r * outNodes.chunk, (r + 1) * outNodes.chunk, inNodes.coordinates.data() + coffset, sBuffer);
		coffset += sBuffer[size + 3];
		sBuffer[size + 4] = insert(offset[1] + foffset, inFaces.ftype.size() - foffset, fdistribution[r], fdistribution[r + 1], inFaces.ftype.data() + foffset, sBuffer);
		sBuffer[size + 5] = edist[foffset + sBuffer[size + 4] + 1] - edist[foffset];
		sBuffer.insert(sBuffer.end(), inFaces.fnodes.begin() + edist[foffset], inFaces.fnodes.begin() + edist[foffset + sBuffer[size + 4] + 1]);
		foffset += sBuffer[size + 4];
		sBuffer[size + 6] = insert(offset[2] + ooffset, inFaces.owner.size() - ooffset, fdistribution[r], fdistribution[r + 1], inFaces.owner.data() + ooffset, sBuffer);
		ooffset += sBuffer[size + 6];
		sBuffer[size + 7] = insert(offset[3] + noffset, inFaces.neighbor.size() - noffset, fdistribution[r], fdistribution[r + 1], inFaces.neighbor.data() + noffset, sBuffer);
		noffset += sBuffer[size + 7];
		sBuffer[size] = sBuffer.size() - size;
	}
	utils::clearVectors(inNodes.coordinates, inFaces.ftype, inFaces.fnodes, inFaces.owner, inFaces.neighbor, edist);
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot balance owner and neighbor data.\n");
	}
	utils::clearVector(sBuffer);

	inFaces.ftype.reserve(fdistribution[info::mpi::rank + 1] - fdistribution[info::mpi::rank]);
	inFaces.fnodes.reserve(5 * (fdistribution[info::mpi::rank + 1] - fdistribution[info::mpi::rank]));
	inFaces.owner.reserve(fdistribution[info::mpi::rank + 1] - fdistribution[info::mpi::rank]);
	inFaces.neighbor.reserve(fdistribution[info::mpi::rank + 1] - fdistribution[info::mpi::rank]);
	outNodes.coordinates.reserve(outNodes.chunk);

	ivector<esint> rbegin(info::mpi::size);
	for (esint r = 0, roffset = 0; r < info::mpi::size; ++r) {
		rbegin[rBuffer[roffset + 2]] = roffset;
		roffset += rBuffer[roffset];
	}
	for (esint r = 0; r < info::mpi::size; ++r) {
		esint offset = rbegin[r];
		esint* data = rBuffer.data() + offset + 8;

		outNodes.coordinates.insert(outNodes.coordinates.end(), reinterpret_cast<_Point<esfloat>*>(data), reinterpret_cast<_Point<esfloat>*>(data) + rBuffer[offset + 3]);
		data += utils::reinterpret_size<esint, _Point<esfloat> >(rBuffer[offset + 3]);
		inFaces.ftype.insert(inFaces.ftype.end(), reinterpret_cast<Element::CODE*>(data), reinterpret_cast<Element::CODE*>(data) + rBuffer[offset + 4]);
		data += utils::reinterpret_size<esint, Element::CODE >(rBuffer[offset + 4]);
		inFaces.fnodes.insert(inFaces.fnodes.end(), data, data + rBuffer[offset + 5]);
		data += rBuffer[offset + 5];
		inFaces.owner.insert(inFaces.owner.end(), data, data + rBuffer[offset + 6]);
		data += rBuffer[offset + 6];
		inFaces.neighbor.insert(inFaces.neighbor.end(), data, data + rBuffer[offset + 7]);
		data += rBuffer[offset + 7];
	}
	utils::clearVector(sBuffer);

	// 2. group faces in order to be able to built elements
	ivector<esint> opermutation(inFaces.owner.size()), npermutation(inFaces.neighbor.size());
	std::iota(opermutation.begin(), opermutation.end(), 0);
	std::sort(opermutation.begin(), opermutation.end(), [&] (const esint &i, const esint &j) { return inFaces.owner[i] != inFaces.owner[j] ? inFaces.owner[i] < inFaces.owner[j] : i < j; });
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (const esint &i, const esint &j) { return inFaces.neighbor[i] != inFaces.neighbor[j] ? inFaces.neighbor[i] < inFaces.neighbor[j] : i < j; });

	edist.reserve(inFaces.ftype.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < inFaces.ftype.size(); ++e) {
		edist.push_back(edist.back() + Element::encode(inFaces.ftype[e]).nodes);
	}

	ivector<esint>::const_iterator oit = opermutation.begin();
	ivector<esint>::const_iterator nit = npermutation.begin();
	sBuffer.reserve(2 * info::mpi::size + 3 * inFaces.ftype.size() + inFaces.fnodes.size());
	std::vector<int> sent(inFaces.ftype.size(), -1);
	for (esint r = 0; r < info::mpi::size; ++r) {
		size_t size = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target

		while (oit != opermutation.end() && inFaces.owner[*oit] < (r + 1) * outFaces.elements.chunk) {
			sent[*oit] = r;
			sBuffer.push_back(edist[*oit + 1] - edist[*oit]);
			sBuffer.insert(sBuffer.end(), inFaces.fnodes.begin() + edist[*oit], inFaces.fnodes.begin() + edist[*oit + 1]);
			sBuffer.push_back(inFaces.owner[*oit]);
			sBuffer.push_back(*oit < (esint)inFaces.neighbor.size() ? inFaces.neighbor[*oit] : -1);
			++oit;
		}
		while (nit != npermutation.end() && inFaces.neighbor[*nit] < (r + 1) * outFaces.elements.chunk) {
			if (sent[*nit] < r) {
				sBuffer.push_back(edist[*nit + 1] - edist[*nit]);
				sBuffer.insert(sBuffer.end(), inFaces.fnodes.begin() + edist[*nit], inFaces.fnodes.begin() + edist[*nit + 1]);
				sBuffer.push_back(inFaces.owner[*nit]);
				sBuffer.push_back(inFaces.neighbor[*nit]);
			}
			++nit;
		}
		sBuffer[size] = sBuffer.size() - size;
	}

	utils::clearVectors(inFaces.ftype, inFaces.fnodes, inFaces.owner, inFaces.neighbor, edist);
	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot exchange faces.\n");
	}
	utils::clearVector(sBuffer);

	// just estimations
	outFaces.ftype.reserve(7 * outFaces.elements.chunk);
	outFaces.fnodes.reserve(5 * 7 * outFaces.elements.chunk);
	outFaces.owner.reserve(7 * outFaces.elements.chunk);
	outFaces.neighbor.reserve(7 * outFaces.elements.chunk);
	for (esint r = 0, offset = 0; r < info::mpi::size; ++r) {
		esint size = offset + rBuffer[offset];
		offset += 2;
		while (offset < size) {
			switch (rBuffer[offset]) {
			case 3: outFaces.ftype.push_back(Element::CODE::TRIANGLE3); break;
			case 4: outFaces.ftype.push_back(Element::CODE::SQUARE4); break;
			default: outFaces.ftype.push_back(Element::CODE::POLYGON);
			}
			outFaces.fnodes.insert(outFaces.fnodes.end(), rBuffer.begin() + offset + 1, rBuffer.begin() + offset + 1 + rBuffer[offset]);
			offset += rBuffer[offset] + 1;
			outFaces.owner.push_back(rBuffer[offset++]);
			outFaces.neighbor.push_back(rBuffer[offset++]);
		}
	}
}

}
}
