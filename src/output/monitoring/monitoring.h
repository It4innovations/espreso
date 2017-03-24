
#ifndef SRC_OUTPUT_MONITORING_MONITORING_H_
#define SRC_OUTPUT_MONITORING_MONITORING_H_

#include <fstream>

#include "../store.h"

namespace espreso {

class Mesh;
class Region;
enum class Property;

namespace output {

enum class Operation {
	AVERAGE,
	MIN,
	MAX
};

struct Monitor {
	size_t printSize;
	Region* region;
	Property property;
	Operation operation;
};

class Monitoring: public Store {

public:
	Monitoring(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);

	void storeSettings(const Step &step) {};
	void storeSettings(size_t steps) {};
	void storeSettings(const std::vector<size_t> &steps) {};

	void storeFETIData(const Step &step, const Instance &instance) {};

	void storeSolution(const Step &step, const std::vector<Solution*> &solution);
	void finalize();

protected:
	const Mesh *_mesh;
	std::ofstream _os;

	std::vector<Monitor> _monitors;
};

}
}


#endif /* SRC_OUTPUT_MONITORING_MONITORING_H_ */
