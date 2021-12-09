
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

namespace espreso {
namespace builder {

void balance(OrderedMeshDatabase &database, OrderedMesh &mesh)
{
	esint total[2] = { (esint)database.coordinates.size(), (esint)database.etype.size() };
	Communication::allReduce(total, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);

	mesh.ntotal = total[0];
	mesh.etotal = total[1];
	mesh.nsize = mesh.nchunk = total[0] / info::mpi::size + ((total[0] % info::mpi::size) ? 1 : 0);
	mesh.esize = mesh.echunk = total[1] / info::mpi::size + ((total[1] % info::mpi::size) ? 1 : 0);
	mesh.noffset = std::min(mesh.nchunk * info::mpi::rank, mesh.ntotal);
	mesh.eoffset = std::min(mesh.echunk * info::mpi::rank, mesh.etotal);
	if (mesh.ntotal <= mesh.noffset + mesh.nsize) {
		mesh.nsize = mesh.ntotal - mesh.noffset;
	}
	if (mesh.etotal <= mesh.eoffset + mesh.esize) {
		mesh.esize = mesh.etotal - mesh.eoffset;
	}

	std::vector<esint> sBuffer, rBuffer, edist;
	edist.reserve(database.etype.size() + 1);
	edist.push_back(0);
	for (size_t e = 0; e < database.etype.size(); ++e) {
		edist.push_back(edist.back() + Mesh::edata[(int)database.etype[e]].nodes);
	}
	std::sort(database.noffsets.begin(), database.noffsets.end());
	std::sort(database.eoffsets.begin(), database.eoffsets.end());

	{ // compute size of the send buffer
		size_t ssize = 0;
		std::vector<OrderedMeshDatabase::Offset>::const_iterator nit = database.noffsets.cbegin();
		std::vector<OrderedMeshDatabase::Offset>::const_iterator eit = database.eoffsets.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			ssize += 4;
			while (OrderedMeshDatabase::chunk(mesh.nchunk, r, database.noffsets, nit, nbegin, nend)) {
				ssize += 2;
				ssize += utils::reinterpret_size<esint, _Point<esfloat> >(nend - nbegin);
			}
			while (OrderedMeshDatabase::chunk(mesh.echunk, r, database.eoffsets, eit, ebegin, eend)) {
				ssize += 3;
				ssize += utils::reinterpret_size<esint, char>(eend - ebegin);
				ssize += edist[eend] - edist[ebegin];
			}
		}
		sBuffer.reserve(ssize);
	}

	{ // build the send buffer
		std::vector<OrderedMeshDatabase::Offset>::const_iterator nit = database.noffsets.cbegin();
		std::vector<OrderedMeshDatabase::Offset>::const_iterator eit = database.eoffsets.cbegin();
		esint nbegin = 0, nend = 0, ebegin = 0, eend = 0;
		for (int r = 0; r < info::mpi::size; ++r) {
			size_t prevsize = sBuffer.size();
			sBuffer.push_back(0); // total size
			sBuffer.push_back(r); // target
			sBuffer.push_back(0); // nodes
			sBuffer.push_back(0); // elements

			while (OrderedMeshDatabase::chunk(mesh.nchunk, r, database.noffsets, nit, nbegin, nend)) {
				sBuffer.push_back(nit->offset + nbegin - nit->start);
				sBuffer.push_back(nend - nbegin);
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(database.coordinates.data() + nbegin), utils::reinterpret_end<esint>(database.coordinates.data() + nbegin, nend - nbegin));
				++sBuffer[prevsize + 2];
			}
			while (OrderedMeshDatabase::chunk(mesh.echunk, r, database.eoffsets, eit, ebegin, eend)) {
				sBuffer.push_back(eit->offset + ebegin - eit->start);
				sBuffer.push_back(eend - ebegin);
				sBuffer.push_back(edist[eend] - edist[ebegin]);
				// it can cause Invalid read message with Valgrind (we touch after etype array but we never use that data)
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(database.etype.data() + ebegin), utils::reinterpret_end<esint>(database.etype.data() + ebegin, eend - ebegin));
				sBuffer.insert(sBuffer.end(), database.enodes.data() + edist[ebegin], database.enodes.data() + edist[eend]);
				++sBuffer[prevsize + 3];
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}

		std::vector<esint>().swap(edist);
		database.clearNodes();
		database.clearElements();
		database.clearValues();
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::internalFailure("Cannot balance parsed data.\n");
	}
	std::vector<esint>().swap(sBuffer);

	mesh.coordinates.resize(mesh.nsize);
	mesh.etype.resize(mesh.esize);

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
			memcpy(mesh.coordinates.data() + coffset - mesh.noffset, reinterpret_cast<_Point<esfloat>*>(rBuffer.data() + offset), csize * sizeof(_Point<esfloat>));
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
	mesh.edist.reserve(mesh.esize + 1);
	mesh.enodes.reserve(enodesTotal);
	std::sort(eintervals.begin(), eintervals.end());
	mesh.edist.push_back(0);
	mesh.dimension = 0;
	for (size_t i = 0; i < eintervals.size(); ++i) {
		esint eoffset = eintervals[i].first;
		esint esize = rBuffer[eintervals[i].second++];
		eintervals[i].second++; // enodes
		memcpy(mesh.etype.data() + eoffset - mesh.eoffset, reinterpret_cast<char*>(rBuffer.data() + eintervals[i].second), esize);
		eintervals[i].second += utils::reinterpret_size<esint, char>(esize);
		for (esint en = 0; en < esize; ++en) {
			for (int nn = 0; nn < Mesh::edata[(int)mesh.etype[eoffset - mesh.eoffset + en]].nodes; ++nn) {
				mesh.enodes.push_back(rBuffer[eintervals[i].second++]);
			}
			mesh.edist.push_back(mesh.enodes.size());
			switch (Mesh::edata[(int)mesh.etype[eoffset - mesh.eoffset + en]].type) {
			case Element::TYPE::POINT:  mesh.dimension = std::max(0, mesh.dimension); break;
			case Element::TYPE::LINE:   mesh.dimension = std::max(1, mesh.dimension); break;
			case Element::TYPE::PLANE:  mesh.dimension = std::max(2, mesh.dimension); break;
			case Element::TYPE::VOLUME: mesh.dimension = std::max(3, mesh.dimension); break;
			}
		}
	}
}


}
}
