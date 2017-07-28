
#ifndef SRC_OUTPUT_STORE_H_
#define SRC_OUTPUT_STORE_H_

#include <cstddef>
#include <vector>

namespace espreso {

struct Step;
struct Solution;
struct Instance;
struct OutputConfiguration;
enum class ElementType;
enum class Property;

namespace output {

class Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	virtual void storeSettings(size_t steps) =0;
	virtual void storeFETIData(const Step &step, const Instance &instance) =0;
	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties) =0;

	virtual void finalize() =0;

	virtual ~Store() {};

protected:
	Store(const OutputConfiguration &configuration): _configuration(configuration) {};

	const OutputConfiguration &_configuration;
};

}
}

#endif /* SRC_OUTPUT_STORE_H_ */
