
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "wrappers/mpi/communication.h"

#include <numeric>

namespace espreso {
namespace builder {

void initializeSequentialFEM(const InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, const TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, int &dimension)
{
	if (info::mpi::size != 1) {
		eslog::internalFailure("usage of a sequential method during a parallel run.\n");
	}

	eslog::info(" == TOTAL NUMBER OF NODES %65d == \n", input.nodes->coordinates.size());
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", input.elements->etype.size());

	clustered.nodes->coordinates.swap(input.nodes->coordinates);
	clustered.elements->etype.swap(input.elements->etype);
	clustered.elements->enodes.swap(input.elements->enodes);

	clustered.elements->edist.reserve(clustered.elements->etype.size() + 1);
	clustered.elements->edist.push_back(0);
	dimension = 0;
	for (size_t e = 0; e < clustered.elements->etype.size(); ++e) {
		if (clustered.elements->etype[e] != Element::CODE::POLYGON && clustered.elements->etype[e] != Element::CODE::POLYHEDRON) {
			clustered.elements->edist.push_back(clustered.elements->edist.back() + Mesh::element(clustered.elements->etype[e]).nodes);
		} else {
			clustered.elements->edist.push_back(clustered.elements->edist.back() + clustered.elements->enodes[clustered.elements->edist.back()]);
		}
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

static void _initializeSequentialFVM(const InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, const TemporalSequentialMesh<ClusteredNodes, OrderedFacesBalanced> &clustered)
{
	if (info::mpi::size != 1) {
		eslog::internalFailure("usage of a sequential method during a parallel run.\n");
	}

	eslog::info(" == TOTAL NUMBER OF NODES %65lu == \n", input.nodes->coordinates.size());
	eslog::info(" == TOTAL NUMBER OF ELEMENTS %62d == \n", input.elements->elements);

	clustered.elements->offset = 0;
	clustered.elements->size = input.elements->elements;
	clustered.elements->total = input.elements->elements;

	clustered.nodes->coordinates.swap(input.nodes->coordinates);
	clustered.elements->etype.swap(input.elements->etype);
	clustered.elements->enodes.swap(input.elements->enodes);
	clustered.elements->owner.swap(input.elements->owner);
	clustered.elements->neighbor.swap(input.elements->neighbor);

	clustered.elements->edist.reserve(clustered.elements->etype.size() + 1);
	clustered.elements->edist.push_back(0);
	for (size_t e = 0; e < clustered.elements->etype.size(); ++e) {
		if (clustered.elements->etype[e] != Element::CODE::POLYGON && clustered.elements->etype[e] != Element::CODE::POLYHEDRON) {
			clustered.elements->edist.push_back(clustered.elements->edist.back() + Mesh::element(clustered.elements->etype[e]).nodes);
		} else {
			clustered.elements->edist.push_back(clustered.elements->edist.back() + clustered.elements->enodes[clustered.elements->edist.back()]);
		}
	}
}

void initializeSequentialFVM(const InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, const TemporalSequentialMesh<ClusteredNodes, OrderedFacesBalanced> &clustered)
{
	_initializeSequentialFVM(input, clustered);

	clustered.nodes->offsets.resize(clustered.nodes->coordinates.size());
	std::iota(clustered.nodes->offsets.begin(), clustered.nodes->offsets.end(), 0);
}

void initializeSequentialFVM(const InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, const TemporalSequentialMesh<ClusteredNodes, OrderedFacesBalanced> &clustered)
{
	InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> _input;
	_input.nodes = input.nodes;
	_input.elements = input.elements;
	_input.regions = input.regions;
	_initializeSequentialFVM(_input, clustered);
	_input.nodes = nullptr; _input.elements = nullptr; _input.regions = nullptr;

	clustered.nodes->offsets.resize(clustered.nodes->coordinates.size());
	for (size_t i = 0; i < input.nodes->offsets.size(); ++i) {
		std::iota(clustered.nodes->offsets.begin() + input.nodes->offsets[i].local, clustered.nodes->offsets.begin() + input.nodes->offsets[i].local + input.nodes->offsets[i].size, input.nodes->offsets[i].global);
	}
}

}
}
