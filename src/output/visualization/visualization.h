
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_

#include "basis/containers/point.h"
#include "output/output.h"

namespace espreso {

struct NamedData;

class Visualization: public OutputWriter {

public:
	static bool isRoot();
	static bool storeData(const NamedData *data, const step::Step &step);
	static Point shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio);

	bool storeStep(const step::Step &step);

	Visualization();
	~Visualization();

};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VISUALIZATION_H_ */
