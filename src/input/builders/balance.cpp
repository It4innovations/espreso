
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
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

void balanceFEM(const InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, int &dimension)
{
	esint total[2] = { (esint)input.nodes->coordinates.size(), (esint)input.elements->etype.size() };
	Communication::allReduce(total, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);

	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", total[0]);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", total[1]);

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
		utils::clearVectors(input.nodes->coordinates, input.nodes->offsets);
		utils::clearVectors(input.elements->etype, input.elements->enodes, input.elements->offsets);
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

void balanceFVM(const InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, const TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> &ordered)
{
	std::vector<esint> sum(4), offset = { (esint)input.nodes->coordinates.size(), (esint)input.elements->etype.size(), (esint)input.elements->owner.size(), (esint)input.elements->neighbor.size() };
	std::vector<esint> fdistribution = Communication::getDistribution<esint>(input.elements->etype.size());
	Communication::exscan(sum, offset);

	distribute(*ordered.nodes, sum[0]);
	distribute(*ordered.elements, input.elements->elements);

	{ // exchange owner, neighbour
		esint osize = input.elements->owner.size(), ostart = offset[2], oend = offset[2] + osize, obegin = 0, olBegin = 0, olEnd = 0;
		esint nsize = input.elements->neighbor.size(), nstart = offset[3], nend = offset[3] + nsize, nbegin = 0, nlBegin = 0, nlEnd = 0;
		esint bsize = 5 * info::mpi::size;
		if (ostart < fdistribution[info::mpi::rank]) {
			bsize += std::min(fdistribution[info::mpi::rank] - ostart, osize);
		}
		if (fdistribution[info::mpi::rank + 1] <= oend) {
			bsize += std::min(oend - fdistribution[info::mpi::rank + 1], osize);
		}
		if (nstart < fdistribution[info::mpi::rank]) {
			bsize += std::min(fdistribution[info::mpi::rank] - nstart, nsize);
		}
		if (fdistribution[info::mpi::rank + 1] <= nend) {
			bsize += std::min(nend - fdistribution[info::mpi::rank + 1], nsize);
		}

		ivector<esint> sBuffer, rBuffer;
		sBuffer.reserve(bsize);
		for (int r = 0; r < info::mpi::size; ++r) {
			size_t size = sBuffer.size();
			sBuffer.push_back(0); // total size
			sBuffer.push_back(r); // rank
			sBuffer.push_back(info::mpi::rank); // rank
			sBuffer.push_back(0); // owners
			sBuffer.push_back(0); // neighbors
			if (r != info::mpi::rank) {
				if (ostart + obegin < fdistribution[r + 1]) {
					sBuffer[size + 3] = std::min(osize, fdistribution[r + 1] - fdistribution[r]);
					sBuffer.insert(sBuffer.end(), input.elements->owner.begin() + obegin, input.elements->owner.begin() + obegin + sBuffer[size + 3]);
					obegin += sBuffer[size + 3]; osize -= sBuffer[size + 3];
				}
				if (nstart + nbegin < fdistribution[r + 1]) {
					sBuffer[size + 4] = std::min(nsize, fdistribution[r + 1] - fdistribution[r]);
					sBuffer.insert(sBuffer.end(), input.elements->neighbor.begin() + nbegin, input.elements->neighbor.begin() + nbegin + sBuffer[size + 4]);
					nbegin += sBuffer[size + 4]; nsize -= sBuffer[size + 4];
				}
			} else {
				olBegin = obegin;
				nlBegin = nbegin;
				obegin += std::min(osize, fdistribution[r + 1] - fdistribution[r]);
				nbegin += std::min(nsize, fdistribution[r + 1] - fdistribution[r]);
				olEnd = obegin;
				nlEnd = nbegin;
				osize -= obegin - olBegin;
				nsize -= nbegin - nlBegin;
			}
			sBuffer[size] = sBuffer.size() - size;
		}
		if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
			eslog::internalFailure("Cannot balance owner and neighbor data.\n");
		}
		utils::clearVector(sBuffer);

		ivector<esint> owner, neighbor;
		owner.reserve(input.elements->etype.size());
		neighbor.reserve(input.elements->etype.size());
		ivector<esint> rbegin(info::mpi::size);
		esint roffset = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			esint start = roffset;
			esint size = rBuffer[roffset++];
			++roffset; // me
			rbegin[rBuffer[roffset++]] = start;
			roffset = start + size;
		}

		for (int r = 0; r < info::mpi::size; ++r) {
			roffset = rbegin[r];
			++roffset; // size
			++roffset; // me
			++roffset; // from
			esint owners = rBuffer[roffset++];
			esint neighbors = rBuffer[roffset++];
			if (r != info::mpi::rank) {
				owner.insert(owner.end(), rBuffer.begin() + roffset, rBuffer.begin() + roffset + owners);
				neighbor.insert(neighbor.end(), rBuffer.begin() + roffset + owners, rBuffer.begin() + roffset + owners + neighbors);
			} else {
				owner.insert(owner.end(), input.elements->owner.begin() + olBegin, input.elements->owner.begin() + olEnd);
				neighbor.insert(neighbor.end(), input.elements->neighbor.begin() + nlBegin, input.elements->neighbor.begin() + nlEnd);
			}
		}
		neighbor.resize(owner.size(), info::mpi::size * ordered.elements->chunk);
		input.elements->owner.swap(owner);
		input.elements->neighbor.swap(neighbor);
	}

	ivector<esint> opermutation(input.elements->owner.size()), npermutation(input.elements->neighbor.size());
	std::iota(opermutation.begin(), opermutation.end(), 0);
	std::sort(opermutation.begin(), opermutation.end(), [&] (const esint &i, const esint &j) { return input.elements->owner[i] != input.elements->owner[j] ? input.elements->owner[i] < input.elements->owner[j] : i < j; });
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (const esint &i, const esint &j) { return input.elements->neighbor[i] != input.elements->neighbor[j] ? input.elements->neighbor[i] < input.elements->neighbor[j] : i < j; });

	ivector<esint> sBuffer, rBuffer, edist;
	edist.reserve(input.elements->etype.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < input.elements->etype.size(); ++e) {
		edist.push_back(edist.back() + Mesh::element(input.elements->etype[e]).nodes);
	}

	{ // compute size of the send buffer
		ivector<esint>::const_iterator oit = opermutation.begin();
		ivector<esint>::const_iterator nit = npermutation.begin();
		std::vector<int> sendto(input.elements->owner.size(), -1);
		size_t ssize = 0;
		if (offset[0] < ordered.nodes->offset) {
			ssize += utils::reinterpret_size<esint, _Point<esfloat> >(std::min(ordered.nodes->offset - offset[0], (esint)input.nodes->coordinates.size()));
		}
		if (ordered.nodes->offset + ordered.nodes->size <= (esint)input.nodes->coordinates.size()) {
			ssize += utils::reinterpret_size<esint, _Point<esfloat> >(std::min(offset[0] + (esint)input.nodes->coordinates.size() - (ordered.nodes->offset + ordered.nodes->size), (esint)input.nodes->coordinates.size()));
		}
		for (int r = 0; r < info::mpi::size; ++r) {
			ssize += 6;
			if (r != info::mpi::rank) {
				while (oit < opermutation.end() && input.elements->owner[*oit] < (r + 1) * ordered.elements->chunk) {
					ssize += 4 + edist[*oit + 1] - edist[*oit]; // owner, neighbor, offset, nodes
					sendto[*oit] = r;
					++oit;
				}
				while (nit < npermutation.end() && input.elements->neighbor[*nit] < (r + 1) * ordered.elements->chunk) {
					if (sendto[*nit] < r) {
						ssize += 4 + edist[*nit + 1] - edist[*nit]; // owner, neighbor, offset, nodes
					}
					++nit;
				}
			} else {
				while (oit < opermutation.end() && input.elements->owner[*oit] < (r + 1) * ordered.elements->chunk) { ++oit; }
				while (nit < npermutation.end() && input.elements->neighbor[*nit] < (r + 1) * ordered.elements->chunk) { ++nit; }
			}
		}
		sBuffer.reserve(ssize);
	}

	esint clBegin = 0, clEnd = 0, oBegin = 0, oEnd = 0, nBegin = 0, nEnd = 0;
	{ // build the send buffer
		ivector<esint>::const_iterator oit = opermutation.begin();
		ivector<esint>::const_iterator nit = npermutation.begin();
		std::vector<int> sendto(input.elements->owner.size(), -1);
		esint csize = input.nodes->coordinates.size(), cstart = offset[0], cbegin = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0); // total size
			sBuffer.push_back(r); // target
			sBuffer.push_back(info::mpi::rank); // me
			sBuffer.push_back(0); // coordinates
			sBuffer.push_back(0); // faces
			sBuffer.push_back(0); // nodes

			if (r != info::mpi::rank) {
				if (cstart + cbegin < (r + 1) * ordered.nodes->chunk) {
					sBuffer[prevsize + 3] = std::min(csize, ordered.nodes->chunk);
					sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(input.nodes->coordinates.data() + cbegin), utils::reinterpret_end<esint>(input.nodes->coordinates.data() + cbegin, sBuffer[prevsize + 3]));
					cbegin += sBuffer[prevsize + 3]; csize -= sBuffer[prevsize + 3];
				}
				auto _oit = oit, _nit = nit;
				while (oit < opermutation.end() && input.elements->owner[*oit] < (r + 1) * ordered.elements->chunk) {
					sBuffer[prevsize + 4] += 1;
					sBuffer[prevsize + 5] += edist[*oit + 1] - edist[*oit];
					sendto[*oit] = r;
					++oit;
				}
				while (nit < npermutation.end() && input.elements->neighbor[*nit] < (r + 1) * ordered.elements->chunk) {
					if (sendto[*nit] < r) {
						sBuffer[prevsize + 4] += 1;
						sBuffer[prevsize + 5] += edist[*nit + 1] - edist[*nit];
					}
					++nit;
				}
				oit = _oit; nit = _nit;
				esint faces = 0, nodes = 0, soffset = sBuffer.size();
				sBuffer.resize(sBuffer.size() + 4 * sBuffer[prevsize + 4] + sBuffer[prevsize + 5]);
				while (oit < opermutation.end() && input.elements->owner[*oit] < (r + 1) * ordered.elements->chunk) {
					sBuffer[soffset + 0 * sBuffer[prevsize + 4] + faces] = input.elements->owner[*oit];
					sBuffer[soffset + 1 * sBuffer[prevsize + 4] + faces] = input.elements->neighbor[*oit];
					sBuffer[soffset + 2 * sBuffer[prevsize + 4] + faces] = offset[1] + *oit;
					sBuffer[soffset + 3 * sBuffer[prevsize + 4] + faces] = edist[*oit + 1] - edist[*oit];
					++faces;
					for (esint n = edist[*oit]; n < edist[*oit + 1]; ++n, ++nodes) {
						sBuffer[soffset + 4 * sBuffer[prevsize + 4] + nodes] = input.elements->enodes[n];
					}
					++oit;
				}
				while (nit < npermutation.end() && input.elements->neighbor[*nit] < (r + 1) * ordered.elements->chunk) {
					if (sendto[*nit] < r) {
						sBuffer[soffset + 0 * sBuffer[prevsize + 4] + faces] = input.elements->owner[*nit];
						sBuffer[soffset + 1 * sBuffer[prevsize + 4] + faces] = input.elements->neighbor[*nit];
						sBuffer[soffset + 2 * sBuffer[prevsize + 4] + faces] = offset[1] + *nit;
						sBuffer[soffset + 3 * sBuffer[prevsize + 4] + faces] = edist[*nit + 1] - edist[*nit];
						++faces;
						for (esint n = edist[*nit]; n < edist[*nit + 1]; ++n, ++nodes) {
							sBuffer[soffset + 4 * sBuffer[prevsize + 4] + nodes] = input.elements->enodes[n];
						}
					}
					++nit;
				}
			} else {
				clBegin = cbegin;
				cbegin += std::min(csize, ordered.nodes->chunk);
				clEnd = cbegin;
				sBuffer[prevsize + 3] = clEnd - clBegin;
				csize -= cbegin - clBegin;
				oBegin = oit - opermutation.begin();
				nBegin = nit - npermutation.begin();
				while (oit < opermutation.end() && input.elements->owner[*oit] < (r + 1) * ordered.elements->chunk) {
					sBuffer[prevsize + 4] += 1;
					sBuffer[prevsize + 5] += edist[*oit + 1] - edist[*oit];
					sendto[*oit] = r;
					++oit;
				}
				while (nit < npermutation.end() && input.elements->neighbor[*nit] < (r + 1) * ordered.elements->chunk) {
					if (sendto[*nit] < r) {
						sBuffer[prevsize + 4] += 1;
						sBuffer[prevsize + 5] += edist[*nit + 1] - edist[*nit];
					}
					++nit;
				}
				oEnd = oit - opermutation.begin();
				nEnd = nit - npermutation.begin();
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot balance parsed data.\n");
	}
	utils::clearVector(sBuffer);

	esint coordinates = 0, faces = 0, nodes = 0, roffset = 0;
	std::vector<esint> rbegin(info::mpi::size); // element offset, rBuffer offset
	for (int r = 0; r < info::mpi::size; ++r) {
		esint start = roffset;
		esint size = rBuffer[roffset++];
		++roffset; // me
		rbegin[rBuffer[roffset++]] = start;
		coordinates += rBuffer[roffset++];
		faces += rBuffer[roffset++];
		nodes += rBuffer[roffset++];
		roffset = start + size;
	}
	ordered.nodes->coordinates.reserve(coordinates);
	ordered.elements->etype.reserve(faces);
	ordered.elements->foffset.reserve(faces);
	ordered.elements->owner.reserve(faces);
	ordered.elements->neighbor.reserve(faces);
	ordered.elements->edist.reserve(faces + 1);
	ordered.elements->edist.push_back(0);
	for (int r = 0; r < info::mpi::size; ++r) {
		if (r != info::mpi::rank) {
			roffset = rbegin[r];
			++roffset; // size
			++roffset; // me
			++roffset; // from
			coordinates = rBuffer[roffset++];
			faces = rBuffer[roffset++];
			nodes = rBuffer[roffset++];
			ordered.nodes->coordinates.insert(ordered.nodes->coordinates.end(), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + roffset), reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + roffset) + coordinates);
			roffset += utils::reinterpret_size<esint, _Point<esfloat> >(coordinates);
			ordered.elements->owner.insert(ordered.elements->owner.end(), rBuffer.begin() + roffset, rBuffer.begin() + roffset + faces);
			roffset += faces;
			ordered.elements->neighbor.insert(ordered.elements->neighbor.end(), rBuffer.begin() + roffset, rBuffer.begin() + roffset + faces);
			roffset += faces;
			ordered.elements->foffset.insert(ordered.elements->foffset.end(), rBuffer.begin() + roffset, rBuffer.begin() + roffset + faces);
			roffset += faces;
			for (esint i = 0; i < faces; ++i, ++roffset) {
				switch (rBuffer[roffset]) {
				case 3: ordered.elements->etype.push_back(Element::CODE::TRIANGLE3); break;
				case 4: ordered.elements->etype.push_back(Element::CODE::SQUARE4); break;
				default: ordered.elements->etype.push_back(Element::CODE::NOT_SUPPORTED);
				}
				ordered.elements->edist.push_back(ordered.elements->edist.back() + rBuffer[roffset]);
			}
			ordered.elements->enodes.insert(ordered.elements->enodes.end(), rBuffer.begin() + roffset, rBuffer.begin() + roffset + nodes);
			roffset += nodes;
		} else {
			std::vector<bool> inserted(input.elements->owner.size(), false);
			ordered.nodes->coordinates.insert(ordered.nodes->coordinates.end(), input.nodes->coordinates.begin() + clBegin, input.nodes->coordinates.begin() + clEnd);
			for (auto it = opermutation.begin() + oBegin; it != opermutation.begin() + oEnd; ++it) {
				ordered.elements->owner.push_back(input.elements->owner[*it]);
				ordered.elements->neighbor.push_back(input.elements->neighbor[*it]);
				ordered.elements->foffset.push_back(offset[1] + *it);
				ordered.elements->etype.push_back(input.elements->etype[*it]);
				ordered.elements->edist.push_back(ordered.elements->edist.back() + edist[*it + 1] - edist[*it]);
				for (esint n = edist[*it]; n < edist[*it + 1]; ++n) {
					ordered.elements->enodes.push_back(input.elements->enodes[n]);
				}
				inserted[*it] = true;
			}
			for (auto it = npermutation.begin() + nBegin; it != npermutation.begin() + nEnd; ++it) {
				if (!inserted[*it]) {
					ordered.elements->owner.push_back(input.elements->owner[*it]);
					ordered.elements->neighbor.push_back(input.elements->neighbor[*it]);
					ordered.elements->foffset.push_back(offset[1] + *it);
					ordered.elements->etype.push_back(input.elements->etype[*it]);
					ordered.elements->edist.push_back(ordered.elements->edist.back() + edist[*it + 1] - edist[*it]);
					for (esint n = edist[*it]; n < edist[*it + 1]; ++n) {
						ordered.elements->enodes.push_back(input.elements->enodes[n]);
					}
				}
			}
		}
	}

	utils::clearVectors(input.nodes->coordinates);
	utils::clearVectors(input.elements->etype, input.elements->enodes, input.elements->owner, input.elements->neighbor);
}

}
}
