
#ifndef SRC_OUTPUT_VISUALIZATION_NETGEN_H_
#define SRC_OUTPUT_VISUALIZATION_NETGEN_H_

#include "visualization.h"
#include "writer/netgenwritter.h"

namespace espreso {

class Mesh;

struct Netgen: public Visualization {
	Netgen();
	~Netgen();

	void updateMesh();
	void updateSolution();

protected:
	NetgenASCIIWritter _writer;
};

}



#endif /* SRC_OUTPUT_VISUALIZATION_NETGEN_H_ */
