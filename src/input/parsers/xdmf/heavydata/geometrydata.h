
#ifndef SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GEOMETRYDATA_H_
#define SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GEOMETRYDATA_H_

#include "basis/containers/allocators.h"

#include <vector>
#include <string>

namespace espreso {

struct MeshBuilder;
class XDMFGeometry;
class XDMFDataItem;
class HDF5;

class GeometryData {
public:
	GeometryData(XDMFGeometry *geometry, XDMFDataItem *geometrydata);
	void read(HDF5 &file);

	int dimension;
	std::string name;
	std::vector<esint> distribution;
	std::vector<float, initless_allocator<float> > data;
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GEOMETRYDATA_H_ */
