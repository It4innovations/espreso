
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

void trivialUpdate(OrderedNodes &ordered, OrderedNodesBalanced &balanced)
{
	swap(balanced, ordered);
	balanced.chunk = balanced.coordinates.size();
	balanced.offset = balanced.chunk;
	balanced.size = balanced.chunk;
	balanced.total = Communication::exscan(balanced.offset);
	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", balanced.total);
}

void trivialUpdate(OrderedNodesBalanced &balanced, ClusteredNodes &clustered)
{
	swap(clustered, balanced);
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

void trivialUpdate(OrderedElements &ordered, OrderedElementsBalanced &balanced)
{
	swap(ordered, balanced);
	balanced.chunk = balanced.etype.size();
	balanced.offset = balanced.chunk;
	balanced.size = balanced.chunk;
	balanced.total = Communication::exscan(balanced.offset);
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", balanced.total);
}

void trivialUpdate(OrderedElementsBalanced &balanced, ClusteredElements &clustered)
{
	swap(balanced, clustered);
	clustered.offsets.resize(balanced.size);
	std::iota(clustered.offsets.begin(), clustered.offsets.end(), balanced.offset);
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


