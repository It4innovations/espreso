
#ifndef SRC_OUTPUT_MONITORING_MONITORING_H_
#define SRC_OUTPUT_MONITORING_MONITORING_H_

#include <fstream>

#include "../store.h"

namespace espreso {

class Mesh;
class Region;
enum class Property;
enum StatisticalData: int;

struct Monitor {
	size_t printSize;
	Region* region;
	std::vector<Property> properties;
	StatisticalData statistics;
};

class Monitoring: public Store {

public:
	Monitoring(const OutputConfiguration &output, const Mesh *mesh);

	void storeSettings(size_t steps) {};
	void storeFETIData(const Step &step, const Instance &instance) {};
	void storeSubSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties) {};

	void storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties);
	void finalize();

	static char delimiter;

protected:
	const Mesh *_mesh;
	std::ofstream _os;


	std::vector<Monitor> _monitors;

private:
	std::vector<Property> getProperties(const std::string &name);
};

}


#endif /* SRC_OUTPUT_MONITORING_MONITORING_H_ */
