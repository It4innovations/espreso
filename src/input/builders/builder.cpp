
#include "builder.h"

#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

namespace espreso {
namespace builder {

struct PackedNodes {
	struct Node {
		Node(char *node)
		: offset(*reinterpret_cast<esint*>(node)),
		  coordinate(*reinterpret_cast<_Point<esfloat>*>(node + sizeof(esint)))
		{

		}

		esint &offset;
		_Point<esfloat> &coordinate;
	};

	std::vector<esint, initless_allocator<esint> > distribution;
	std::vector<char, initless_allocator<char> > data;
	size_t nodes, regions = 0, values = 0;

	Node operator[](size_t n) { return Node{data.data() + n * _size}; }

	// ID, coordinates, regions bit-array, value
	void set(size_t nodes, size_t regions, size_t values)
	{
		this->nodes = nodes;
		this->regions = regions;
		this->values = values;
		_size = sizeof(esint) + 3 * sizeof(esfloat) + regions + values;
		data.resize(nodes * _size);
	}

	void print()
	{
		Communication::serialize([&] () {
			printf(" === rank %d ===\n", info::mpi::rank);
			for (size_t n = 0; n < nodes; ++n) {
				Node node{data.data() + n * _size};
				printf("%d: %f %f %f\n", node.offset, node.coordinate.x, node.coordinate.y, node.coordinate.z);
			}
		});
	}

protected:
	size_t _size;
};

struct PackedElements {
	struct Element {
		Element(char *element)
		: offset(*reinterpret_cast<esint*>(element)),
		  type(*reinterpret_cast<int*>(element + sizeof(esint))),
		  nodes(*reinterpret_cast<int*>(element + sizeof(esint) + sizeof(int))),
		  node(reinterpret_cast<esint*>(element + sizeof(esint) + sizeof(int) + sizeof(int)))
		{

		}

		esint &offset;
		int &type;
		int &nodes;
		esint *node;
	};

	std::vector<esint, initless_allocator<esint> > distribution;
	std::vector<char, initless_allocator<char> > data;
	size_t elements, regions, values;

	Element operator[](size_t e) { return Element{data.data() + distribution[e]}; }

	void set(const std::vector<esint> &esize, size_t regions, size_t values)
	{
		this->elements = esize.size();
		this->regions = regions;
		this->values = values;
		distribution.resize(esize.size() + 1);
		distribution[0] = 0;
		for (size_t e = 0; e < esize.size(); ++e) {
			distribution[e + 1] = distribution[e] + sizeof(esint) + sizeof(int) + sizeof(int) + esize[e] * sizeof(esint) + regions + values;
		}
		data.resize(distribution.back());
	}

	void print()
	{
		Communication::serialize([&] () {
			printf(" === rank %d ===\n", info::mpi::rank);
			for (size_t e = 0; e < elements; ++e) {
				Element element{data.data() + distribution[e]};
				printf("%d: %d %d ->", element.offset, element.type, element.nodes);
				for (int n = 0; n < element.nodes; ++n) {
					printf(" %d", element.node[n]);
				}
				printf("\n");
			}
		});
	}
};

struct _Mesh {
	size_t nchunk, nsize;
	size_t echunk, esize;

	PackedNodes nodes;
	PackedElements elements;
};

static void connect(OrderedMeshDatabase &database, _Mesh &mesh)
{
	if (database.noffset.size() != database.coordinates.size()) {
		eslog::internalFailure("the size of OrderedMeshDatabase::noffset != OrderedMeshDatabase::coordinates.\n");
	}

	mesh.nodes.set(database.coordinates.size(), 0, 0);

	#pragma omp parallel for
	for (size_t n = 0; n < database.noffset.size(); ++n) {
		mesh.nodes[n].offset = database.noffset[n];
		mesh.nodes[n].coordinate = database.coordinates[n];
	}
	database.clearNodes();

	if (database.eoffset.size() != database.esize.size()) {
		eslog::internalFailure("the size of OrderedMeshDatabase::eoffset != OrderedMeshDatabase::esize.\n");
	}
	if (database.eoffset.size() != database.etype.size()) {
		eslog::internalFailure("the size of OrderedMeshDatabase::eoffset != OrderedMeshDatabase::etype.\n");
	}
	mesh.elements.set(database.esize, 0, 0);

	for (size_t e = 0, edist = 0; e < database.eoffset.size(); edist += database.esize[e++]) {
		mesh.elements[e].offset = database.eoffset[e];
		mesh.elements[e].type = database.etype[e];
		mesh.elements[e].nodes = database.esize[e];
		for (int n = 0; n < mesh.elements[e].nodes; ++ n) {
			mesh.elements[e].node[n] = database.enodes[edist + n];
		}
	}
	database.clearElements();

	mesh.nodes.print();
	mesh.elements.print();

	Communication::serialize([&] () {
		printf(" === rank %d ===\n", info::mpi::rank);
		printf(" === nvalues ===\n");
		for (size_t i = 0; i < database.nvalues.size(); ++i) {
			for (size_t j = 0; j < database.nvalues[i].values.size(); ++j) {
				for (size_t k = 0; k < database.nvalues[i].values[j].values.size(); ++k) {
					printf("%d: %f\n", database.nvalues[i].values[j].begin + k, database.nvalues[i].values[j].values[k]);
				}
			}
		}
		printf(" === evalues ===\n");
		for (size_t i = 0; i < database.evalues.size(); ++i) {
			for (size_t j = 0; j < database.evalues[i].values.size(); ++j) {
				for (size_t k = 0; k < database.evalues[i].values[j].values.size(); ++k) {
					printf("%d: %f\n", database.evalues[i].values[j].begin + k, database.evalues[i].values[j].values[k]);
				}
			}
		}
	});

	size_t totalsize[2] = { mesh.nodes.nodes, mesh.elements.elements };
	Communication::allReduce(totalsize, NULL, 2, MPITools::getType<size_t>().mpitype, MPI_SUM);

	mesh.nsize = totalsize[0];
	mesh.esize = totalsize[1];
	mesh.nchunk = mesh.nsize / info::mpi::size + (mesh.nsize % info::mpi::size ? 1 : 0);
	mesh.echunk = mesh.esize / info::mpi::size + (mesh.esize % info::mpi::size ? 1 : 0);

//	std::vector<esint> sBuffer, rBuffer;
//	sBuffer.reserve(3 * info::mpi::size + _meshData.nIDs.size() + _meshData.coordinates.size() * sizeof(Point) / sizeof(esint));
//	size_t prevsize;
//	auto nbegin = permutation.begin();
//	for (int r = 0; r < info::mpi::size; r++) {
//		prevsize = sBuffer.size();
//		sBuffer.push_back(0); // total size
//		sBuffer.push_back(r); // target
//		sBuffer.push_back(0); // number of coordinates
//
//		auto n = nbegin;
//		for ( ; n != permutation.end() && _meshData.nIDs[*n] < _nDistribution[r + 1]; ++n) {
//			sBuffer.push_back(_meshData.nIDs[*n]);
//			sBuffer.insert(sBuffer.end(), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n + 1));
//		}
//		sBuffer[prevsize + 2] = n - nbegin;
//		nbegin = n;
//
//		sBuffer[prevsize] = sBuffer.size() - prevsize;
//	}
}

static void buildSequential(OrderedMeshDatabase &database, Mesh &mesh)
{
	_Mesh internal;

	connect(database, internal);
}

void build(OrderedMeshDatabase &database, Mesh &mesh)
{
	if (info::mpi::size == 1) {
		buildSequential(database, mesh);
		return;
	}

	eslog::startln("BUILDER: BUILD SCATTERED MESH", "BUILDER");

	_Mesh internal;

	connect(database, internal);
	eslog::checkpointln("BUILDER: DATA CONNECTED");

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
