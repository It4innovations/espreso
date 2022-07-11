
#include "builder.utils.h"
#include "esinfo/eslog.hpp"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

void swap(Nodes &n1, Nodes &n2)
{
	n1.coordinates.swap(n2.coordinates);
}

void swap(Elements &e1, Elements &e2)
{
	e1.etype.swap(e2.etype);
	e1.enodes.swap(e2.enodes);
}

void swap(Faces &f1, Faces &f2)
{
	f1.ftype.swap(f2.ftype);
	f1.fnodes.swap(f2.fnodes);
	f1.owner.swap(f2.owner);
	f1.neighbor.swap(f2.neighbor);
}

void trivialUpdate(NodesBlocks &blocks, OrderedNodesChunked &chunked)
{
	swap(chunked, blocks);
	chunked.offset = chunked.coordinates.size();
	chunked.size = chunked.offset;
	chunked.total = Communication::exscan(chunked.offset);
	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", chunked.total);
}

void trivialUpdate(FacesBlocks &blocks, OrderedFacesChunked &chunked)
{
	swap(chunked, blocks);
	chunked.elements.offset = 0;
	for (size_t i = 0; i < blocks.elements.blocks.size(); ++i) {
		chunked.elements.offset += blocks.elements.blocks[i].size;
	}
	chunked.elements.size = chunked.elements.offset;
	chunked.elements.total = Communication::exscan(chunked.elements.offset);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", chunked.elements.total);
}

void trivialUpdate(NodesBlocks &blocks, OrderedNodesBalanced &balanced)
{
	swap(balanced, blocks);
	balanced.chunk = balanced.coordinates.size();
	balanced.offset = balanced.chunk;
	balanced.size = balanced.chunk;
	balanced.total = Communication::exscan(balanced.offset);
	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", balanced.total);
}

void trivialUpdate(OrderedNodesChunked &chunked, ClusteredNodes &clustered)
{
	swap(clustered, chunked);
	clustered.offsets.resize(chunked.size);
	std::iota(clustered.offsets.begin(), clustered.offsets.end(), chunked.offset);
}

void trivialUpdate(OrderedNodesBalanced &balanced, ClusteredNodes &clustered)
{
	swap(clustered, balanced);
	clustered.offsets.resize(balanced.size);
	std::iota(clustered.offsets.begin(), clustered.offsets.end(), balanced.offset);
}

void trivialUpdate(OrderedElementsChunked &chunked, ClusteredElements &clustered)
{
	swap(chunked, clustered);
	clustered.offsets.resize(chunked.size);
	std::iota(clustered.offsets.begin(), clustered.offsets.end(), chunked.offset);
}

void trivialUpdate(OrderedElementsBalanced &balanced, ClusteredElements &clustered)
{
	swap(balanced, clustered);
	clustered.offsets.resize(balanced.size);
	std::iota(clustered.offsets.begin(), clustered.offsets.end(), balanced.offset);
}

void trivialUpdate(ClusteredNodes &clustered, MergedNodes &merged)
{
	swap(merged, clustered);
	merged.offsets.swap(clustered.offsets);
	merged.duplication.clear();
}

void trivialUpdate(MergedNodes &merged, LinkedNodes &linked)
{
	swap(linked, merged);
	linked.offsets.swap(merged.offsets);
	linked.duplication.swap(merged.duplication);
	linked.neighbors.clear();
	linked.rankDistribution.resize(linked.offsets.size() + 1);
	std::iota(linked.rankDistribution.begin(), linked.rankDistribution.end(), 0);
	linked.rankData.resize(linked.offsets.size(), info::mpi::rank);
}

void trivialUpdate(ElementsBlocks &blocks, OrderedElementsChunked &chunked)
{
	swap(blocks, chunked);
	chunked.offset = chunked.etype.size();
	chunked.size = chunked.etype.size();
	chunked.total = Communication::exscan(chunked.offset);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", chunked.total);
}

void trivialUpdate(ElementsBlocks &blocks, OrderedElementsBalanced &balanced)
{
	swap(blocks, balanced);
	balanced.chunk = balanced.etype.size();
	balanced.offset = balanced.chunk;
	balanced.size = balanced.chunk;
	balanced.total = Communication::exscan(balanced.offset);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", balanced.total);
}

void trivialUpdate(ClusteredElements &clustered, MergedElements &merged)
{
	swap(clustered, merged);
	swap(clustered, merged);
	clustered.offsets.swap(merged.offsets);
	merged.duplication.clear();
}

}
}


