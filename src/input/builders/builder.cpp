
#include "builder.h"

#include "basis/utilities/packing.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

namespace espreso {
namespace builder {

struct OrderedMesh: public OrderedMeshDatabase {
	esint nchunk, noffset, nsize, ntotal;
	esint echunk, eoffset, esize, etotal;

	static void edist(std::vector<Element::CODE> &ecodes, std::vector<esint> &dist)
	{
		dist.push_back(0);
		for (size_t e = 0; e < ecodes.size(); ++e) {
			dist.push_back(dist.back() + Mesh::edata[(int)ecodes[e]].nodes);
		}
	}
};

static void balance(OrderedMeshDatabase &database, OrderedMesh &mesh)
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
	OrderedMesh::edist(database.etype, edist);
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
				ssize += 2;
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
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(database.etype.data() + ebegin), utils::reinterpret_end<esint>(database.etype.data() + ebegin, eend - ebegin));
				sBuffer.insert(sBuffer.end(), reinterpret_cast<esint*>(database.enodes.data() + edist[ebegin]), utils::reinterpret_end<esint>(database.enodes.data() + edist[ebegin], edist[eend] - edist[ebegin]));
				++sBuffer[prevsize + 3];
			}
			sBuffer[prevsize] = sBuffer.size() - prevsize;
		}

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

	size_t offset = 0;
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
			esint esize = rBuffer[offset++];
			memcpy(mesh.etype.data() + eoffset - mesh.eoffset, reinterpret_cast<char*>(rBuffer.data() + offset), esize);
			offset += utils::reinterpret_size<esint, char>(esize);
			for (esint en = 0; en < esize; ++en) {
				for (int nn = 0; nn < Mesh::edata[(int)mesh.etype[eoffset - mesh.eoffset + en]].nodes; ++nn) {
					mesh.enodes.push_back(rBuffer[offset++]);
				}
			}
		}
	}
}

static void buildSequential(OrderedMeshDatabase &database, Mesh &mesh)
{
	OrderedMesh internal;

	balance(database, internal);
//	connect(database, internal);
}

void build(OrderedMeshDatabase &database, Mesh &mesh)
{
	if (info::mpi::size == 1) {
		buildSequential(database, mesh);
		return;
	}

	eslog::startln("BUILDER: BUILD SCATTERED MESH", "BUILDER");

	OrderedMesh internal;

	balance(database, internal);
	eslog::checkpointln("BUILDER: DATA BALANCED");

	MPI_Finalize();
	exit(0);

////	balance();
////	eslog::checkpointln("BUILDER: DATA BALANCED");
//
//	assignRegions(_meshData.eregions, _meshData.eIDs, _eDistribution, _eregsize, _eregions);
//	assignRegions(_meshData.nregions, _meshData.nIDs, _nDistribution, _nregsize, _nregions);
//	eslog::checkpointln("BUILDER: REGION ASSIGNED");
//
////	reindexRegions();
//
//	assignNBuckets();
//	eslog::checkpointln("BUILDER: NODES BUCKETS COMPUTED");
//
//	assignEBuckets();
//	eslog::checkpointln("BUILDER: ELEMENTS BUCKETS ASSIGNED");
//
//	clusterize();
//	eslog::checkpointln("BUILDER: ELEMENTS CLUSTERED");
//
//	computeSFCNeighbors();
//	eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");
//
//	if (_meshData.removeDuplicates) {
//		mergeDuplicatedNodes();
//		eslog::checkpointln("BUILDER: DUPLICATED NODES MERGED");
//	}
//
//	sortElements();
//	eslog::checkpointln("BUILDER: ELEMENTS SORTED");
//
//	linkup();
//	fillNeighbors();
//	eslog::checkpointln("BUILDER: LINKED UP");
//
//	sortNodes();
//	eslog::checkpointln("BUILDER: NODES SORTED");
//
//	reindexElementNodes();
//	eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");
//
//	if (_meshData.removeDuplicates) {
//		removeDuplicateElements();
//		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS REMOVED");
//	}
//
//	fillNodes();
//	eslog::checkpointln("BUILDER: NODES FILLED");
//
//	fillElements();
//	eslog::checkpointln("BUILDER: ELEMENTS SORTED");
//
//	if (info::mesh->nodes->elements == NULL) {
//		mesh::linkNodesAndElements(info::mesh->elements, info::mesh->nodes, info::mesh->neighbors);
//	}
//
//	exchangeBoundary();
//	eslog::checkpointln("BUILDER: BOUNDARY EXCHANGED");
//
//	fillRegions(_meshData.eregions, _eregsize, _eregions);
//	fillRegions(_meshData.nregions, _nregsize, _nregions);
//	fillElementRegions();
//	fillBoundaryRegions();
//	fillNodeRegions();
//	eslog::checkpointln("BUILDER: REGIONS FILLED");
//
//	reindexBoundaryNodes();
//	eslog::endln("BUILDER: BOUNDARY NODES REINDEXED");
}

}
}
