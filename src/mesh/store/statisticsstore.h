
#ifndef SRC_MESH_STORE_STATISTICSSTORE_H_
#define SRC_MESH_STORE_STATISTICSSTORE_H_

namespace espreso {

struct Statistics {
	double min, max, avg, norm, absmin, absmax;

	Statistics();
	void reset();
};

}

#endif /* SRC_MESH_STORE_STATISTICSSTORE_H_ */
