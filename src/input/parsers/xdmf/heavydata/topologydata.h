
#ifndef SRC_INPUT_FORMATS_XDMF_HEAVYDATA_TOPOLOGYDATA_H_
#define SRC_INPUT_FORMATS_XDMF_HEAVYDATA_TOPOLOGYDATA_H_

#include "basis/containers/allocators.h"

#include <string>
#include <vector>

namespace espreso {

struct MeshBuilder;
class XDMFTopology;
class XDMFDataItem;
class HDF5;

class TopologyData {
public:
	static int firstrank;
	static const int align = 20;

	TopologyData(XDMFTopology *topology, XDMFDataItem *topologydata);
	void read(HDF5 &file);

	int etype, esize;
	std::string name;
	std::vector<esint> distribution;
	std::vector<esint, initless_allocator<esint> > data;
};

}



#endif /* SRC_INPUT_FORMATS_XDMF_HEAVYDATA_TOPOLOGYDATA_H_ */
