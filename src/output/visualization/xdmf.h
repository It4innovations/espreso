
#ifndef SRC_OUTPUT_VISUALIZATION_XDMF_H_
#define SRC_OUTPUT_VISUALIZATION_XDMF_H_

#include "visualization.h"

#include <vector>

namespace espreso {

class HDF5;
class XML;
struct XDMFData;

class XDMF: public Visualization {

public:
	struct HDF5MetaData {
		std::string name;
		int dimension;
		esint offset, size, totalsize;
	};

	struct Geometry: public HDF5MetaData {
		std::vector<float> coordinates;
	};
	struct Topology: public HDF5MetaData {
		std::vector<int> topology;
	};
	struct Attribute: public HDF5MetaData {
		std::vector<float> values;
	};

	XDMF();
	~XDMF();

	void updateMesh();
	void updateMonitors(const step::Step &step);
	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

protected:
	void updateSolution(const step::Step &step);

	HDF5 *_hdf5;
	XML *_xml;
	XDMFData *_data;
};

}




#endif /* SRC_OUTPUT_VISUALIZATION_XDMF_H_ */
