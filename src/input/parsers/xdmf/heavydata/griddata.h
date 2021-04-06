
#ifndef SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GRIDDATA_H_
#define SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GRIDDATA_H_

#include "geometrydata.h"
#include "topologydata.h"

#include <vector>

namespace espreso {

class LightData;
class XDMFGrid;
class XDMFGeometry;
class XDMFTopology;
struct MeshBuilder;

class GridData {
	struct _Grid {
		XDMFGrid *grid;
		XDMFGeometry *geometry;
		XDMFDataItem *geometrydata;
		XDMFTopology *topology;
		XDMFDataItem *topologydata;
	};

public:
	GridData(LightData &root);

	void scan();
	void read();
	void align();
	void parse(MeshBuilder &mesh);

protected:
	LightData &_lightdata;
	std::vector<_Grid> _grids;
	std::vector<GeometryData> _geometry;
	std::vector<TopologyData > _topology;
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_HEAVYDATA_GRIDDATA_H_ */
