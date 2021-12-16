
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

static bool chunk(const esint &mpichunk, const int &rank, const std::vector<DatabaseOffset> &offsets, std::vector<DatabaseOffset>::const_iterator &it, esint &begin, esint &end)
{
	if (it != offsets.end() && end == it->local + it->size) {
		++it;
	}
	if (it == offsets.end()) {
		return false;
	}
	if (it->global + end - it->local == mpichunk * (rank + 1)) {
		return false;
	}
	begin = end = 0;
	while (it != offsets.end() && it->global + it->size < mpichunk * rank) { ++it; }
	if (it != offsets.end() && it->global < mpichunk * (rank + 1)) {
		begin = std::max(it->global, mpichunk * rank);
		end = std::min(it->global + it->size, mpichunk * (rank + 1));
		begin = it->local + begin - it->global;
		end = it->local + end - it->global;
	}
	return begin != end;
}

static void distribute(OrderedDataDistribution &distribution, const esint &total)
{
	distribution.total = total;
	distribution.size = distribution.chunk = total / info::mpi::size + ((total % info::mpi::size) ? 1 : 0);
	distribution.offset = std::min(distribution.chunk * info::mpi::rank, distribution.total);
	if (distribution.total <= distribution.offset + distribution.size) {
		distribution.size = distribution.total - distribution.offset;
	}
}

void balance(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, int &dimension)
{
	esint total[2] = { (esint)input.nodes->coordinates.size(), (esint)input.elements->etype.size() };
	Communication::allReduce(total, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);

	distribute(*ordered.nodes, total[0]);
	distribute(*ordered.elements, total[1]);

	std::vector<esint> sBuffer, rBuffer, edist;
	edist.reserve(input.elements->etype.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < input.elements->etype.size(); ++e) {
		edist.push_back(edist.back() + Mesh::element(input.elements->etype[e]).nodes);
	}
	std::sort(input.nodes->offsets.begin(), input.nodes->offsets.end());
	std::sort(input.elements->offsets.begin(), input.elements->offsets.end());

	{ // compute size of the send buffer
		size_t ssize = 0;
		std::vector<DatabaseOffset>::const_iterator nit = input.nodes->offsets.cbegin();
		std::vector<DatabaseOffset>::const_iterator eit = input.elements->offsets.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			ssize += 4;
			while (chunk(ordered.nodes->chunk, r, input.nodes->offsets, nit, nbegin, nend)) {
				ssize += 2;
				ssize += utils::reinterpret_size<esint, _Point<esfloat> >(nend - nbegin);
			}
			while (chunk(ordered.elements->chunk, r, input.elements->offsets, eit, ebegin, eend)) {
				ssize += 3;
				ssize += utils::reinterpret_size<esint, char>(eend - ebegin);
				ssize += edist[eend] - edist[ebegin];
			}
		}
		sBuffer.reserve(ssize);
	}

	{ // build the send buffer
		std::vector<DatabaseOffset>::const_iterator nit = input.nodes->offsets.cbegin();
		std::vector<DatabaseOffset>::const_iterator eit = input.elements->offsets.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0); // total size
			sBuffer.push_back(r); // target
			sBuffer.push_back(0); // nodes
			sBuffer.push_back(0); // elements

			while (chunk(ordered.nodes->chunk, r, input.nodes->offsets, nit, nbegin, nend)) {
				sBuffer.push_back(nit->global + nbegin - nit->local);
				sBuffer.push_back(nend - nbegin);
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(input.nodes->coordinates.data() + nbegin), utils::reinterpret_end<esint>(input.nodes->coordinates.data() + nbegin, nend - nbegin));
				++sBuffer[prevsize + 2];
			}
			while (chunk(ordered.elements->chunk, r, input.elements->offsets, eit, ebegin, eend)) {
				sBuffer.push_back(eit->global + ebegin - eit->local);
				sBuffer.push_back(eend - ebegin);
				sBuffer.push_back(edist[eend] - edist[ebegin]);
				// it can cause Invalid read message with Valgrind (we touch after etype array but we never use that data)
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(input.elements->etype.data() + ebegin), utils::reinterpret_end<esint>(input.elements->etype.data() + ebegin, eend - ebegin));
				sBuffer.insert(sBuffer.end(), input.elements->enodes.data() + edist[ebegin], input.elements->enodes.data() + edist[eend]);
				++sBuffer[prevsize + 3];
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}

		utils::clearVector(edist);
		delete input.nodes; input.nodes = nullptr;
		delete input.elements; input.elements = nullptr;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot balance parsed data.\n");
	}
	utils::clearVector(sBuffer);

	ordered.nodes->coordinates.resize(ordered.nodes->size);
	ordered.elements->etype.resize(ordered.elements->size);

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
			memcpy(ordered.nodes->coordinates.data() + coffset - ordered.nodes->offset, reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset), csize * sizeof(_Point<esfloat>));
			offset += utils::reinterpret_size<esint, _Point<esfloat> >(csize);
		}
		for (esint e = 0; e < eblocks; ++e) {
			esint eoffset = rBuffer[offset++];
			eintervals.push_back(std::make_pair(eoffset, offset));
			esint esize = rBuffer[offset++];
			esint enodes = rBuffer[offset++];
			offset += utils::reinterpret_size<esint, char>(esize) + enodes;
			enodesTotal += enodes;
		}
	}
	ordered.elements->edist.reserve(ordered.elements->size + 1);
	ordered.elements->enodes.reserve(enodesTotal);
	std::sort(eintervals.begin(), eintervals.end());
	ordered.elements->edist.push_back(0);
	dimension = 0;
	for (size_t i = 0; i < eintervals.size(); ++i) {
		esint eoffset = eintervals[i].first;
		esint esize = rBuffer[eintervals[i].second++];
		eintervals[i].second++; // enodes
		memcpy(ordered.elements->etype.data() + eoffset - ordered.elements->offset, reinterpret_cast<char*>(rBuffer.data() + eintervals[i].second), esize);
		eintervals[i].second += utils::reinterpret_size<esint, char>(esize);
		for (esint en = 0; en < esize; ++en) {
			for (int nn = 0; nn < Mesh::element(ordered.elements->etype[eoffset - ordered.elements->offset + en]).nodes; ++nn) {
				ordered.elements->enodes.push_back(rBuffer[eintervals[i].second++]);
			}
			ordered.elements->edist.push_back(ordered.elements->enodes.size());
			switch (Mesh::element(ordered.elements->etype[eoffset - ordered.elements->offset + en]).type) {
			case Element::TYPE::POINT:  dimension = std::max(0, dimension); break;
			case Element::TYPE::LINE:   dimension = std::max(1, dimension); break;
			case Element::TYPE::PLANE:  dimension = std::max(2, dimension); break;
			case Element::TYPE::VOLUME: dimension = std::max(3, dimension); break;
			}
		}
	}
}


}
}
