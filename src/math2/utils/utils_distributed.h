
#ifndef SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_
#define SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_

#include <vector>

namespace espreso {

struct DataSynchronization {
	std::vector<std::vector<double> > sBuffer, rBuffer;
	std::vector<std::vector<esint> > rIndices;
	std::vector<esint> nStart;

	void gatherFromUpper(double *vals);
	void scatterToUpper(double *vals);
};

}

#endif /* SRC_MATH2_UTILS_UTILS_DISTRIBUTED_H_ */
