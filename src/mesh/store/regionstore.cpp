
#include "regionstore.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "mesh/element.h"

using namespace espreso;

RegionStore::RegionStore(const std::string &name)
: name(name),
  processPerCode(static_cast<int>(Element::CODE::SIZE)), nodes(NULL)
{

}

RegionStore::RegionStore(const char* &packedData)
: RegionStore("")
{
	unpackFull(packedData);
}

RegionStore::~RegionStore()
{
	if (nodes != NULL) { delete nodes; }
}

esint RegionStore::getPosition(esint node) const
{
	return nodeInfo.position[std::lower_bound(nodes->datatarray().begin(), nodes->datatarray().end(), node) - nodes->datatarray().begin()];
}

size_t RegionStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(name);
	packedSize += utils::packedSize(nodes);
	packedSize += utils::packedSize(nodeInfo.nhalo);
	packedSize += utils::packedSize(nodeInfo.offset);
	packedSize += utils::packedSize(nodeInfo.size);
	packedSize += utils::packedSize(nodeInfo.totalSize);
	packedSize += utils::packedSize(nodeInfo.position);

	return packedSize;
}

void RegionStore::packFull(char* &p) const
{
	utils::pack(name, p);
	utils::pack(nodes, p);
	utils::pack(nodeInfo.nhalo, p);
	utils::pack(nodeInfo.offset, p);
	utils::pack(nodeInfo.size, p);
	utils::pack(nodeInfo.totalSize, p);
	utils::pack(nodeInfo.position, p);
}

void RegionStore::unpackFull(const char* &p)
{
	utils::unpack(name, p);
	utils::unpack(nodes, p);
	utils::unpack(nodeInfo.nhalo, p);
	utils::unpack(nodeInfo.offset, p);
	utils::unpack(nodeInfo.size, p);
	utils::unpack(nodeInfo.totalSize, p);
	utils::unpack(nodeInfo.position, p);
}

size_t RegionStore::packedSize() const
{
	return
			utils::packedSize(name) +
			nodes->packedSize() +
			utils::packedSize(nodeInfo.nhalo) +
			utils::packedSize(nodeInfo.offset) +
			utils::packedSize(nodeInfo.size) +
			utils::packedSize(nodeInfo.totalSize) +
			utils::packedSize(nodeInfo.position);
}

void RegionStore::pack(char* &p) const
{
	utils::pack(name, p);
	nodes->pack(p);
	utils::pack(nodeInfo.nhalo, p);
	utils::pack(nodeInfo.offset, p);
	utils::pack(nodeInfo.size, p);
	utils::pack(nodeInfo.totalSize, p);
	utils::pack(nodeInfo.position, p);
}

void RegionStore::unpack(const char* &p)
{
	utils::unpack(name, p);
	if (nodes == NULL) {
		nodes = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	nodes->unpack(p);
	utils::unpack(nodeInfo.nhalo, p);
	utils::unpack(nodeInfo.offset, p);
	utils::unpack(nodeInfo.size, p);
	utils::unpack(nodeInfo.totalSize, p);
	utils::unpack(nodeInfo.position, p);
}


