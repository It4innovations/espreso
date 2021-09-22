
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
	void updateMonitors(step::TYPE type);
	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

protected:
	NetgenASCIIWritter _writer;
};

}



#endif /* SRC_OUTPUT_VISUALIZATION_NETGEN_H_ */
