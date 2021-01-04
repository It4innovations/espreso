
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_

#include "visualization.h"
#include "writer/stlwritter.h"

namespace espreso {

class Mesh;

struct STL: public Visualization {
	STL();
	~STL();

	void updateMesh();
	void updateSolution();

protected:
	STLBinaryOutputWriter _writer;
};

}


#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_COLLECTED_STL_H_ */
