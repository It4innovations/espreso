
#ifndef SRC_OUTPUT_VISUALIZATION_OPENVDB_H_
#define SRC_OUTPUT_VISUALIZATION_OPENVDB_H_

#include "visualization.h"

#include <string>

namespace espreso {

class OpenVDB: public Visualization {
public:
	OpenVDB();
	~OpenVDB();

	void updateMesh();
	void updateSolution();

protected:
	std::string _filename;
};

}

#endif