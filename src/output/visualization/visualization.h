
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_

#include "basis/containers/point.h"
#include "output/resultstore.h"

namespace espreso {

class NamedData;

class Visualization: public ResultStoreBase {

public:
	static bool isRoot();
	static bool storeStep();
	static bool storeData(const NamedData *data);
	static Point shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio);

	void updateMonitors() {}

	Visualization(const Mesh &mesh);
	~Visualization();

};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_ */
