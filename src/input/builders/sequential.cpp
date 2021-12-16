
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

void initialize(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, int &dimension)
{
	if (info::mpi::size != 1) {
		eslog::internalFailure("usage of a sequential method during a parallel run.\n");
	}

	clustered.nodes->coordinates.swap(input.nodes->coordinates);
	clustered.elements->etype.swap(input.elements->etype);
	clustered.elements->enodes.swap(input.elements->enodes);

	clustered.elements->edist.reserve(clustered.elements->etype.size() + 1);
	clustered.elements->edist.push_back(0);
	dimension = 0;
	for (size_t e = 0; e < clustered.elements->etype.size(); ++e) {
		clustered.elements->edist.push_back(clustered.elements->edist.back() + Mesh::element(clustered.elements->etype[e]).nodes);
		switch (Mesh::element(clustered.elements->etype[e]).type) {
		case Element::TYPE::POINT:  dimension = std::max(0, dimension); break;
		case Element::TYPE::LINE:   dimension = std::max(1, dimension); break;
		case Element::TYPE::PLANE:  dimension = std::max(2, dimension); break;
		case Element::TYPE::VOLUME: dimension = std::max(3, dimension); break;
		}
	}

	clustered.nodes->offsets.resize(clustered.nodes->coordinates.size());
	clustered.elements->offsets.resize(clustered.elements->etype.size());
	std::iota(clustered.nodes->offsets.begin(), clustered.nodes->offsets.end(), 0);
	std::iota(clustered.elements->offsets.begin(), clustered.elements->offsets.end(), 0);
}

}
}
