
#ifndef SRC_OUTPUT_STORE_H_
#define SRC_OUTPUT_STORE_H_

#include <vector>

namespace espreso {

struct Step;
struct Solution;
struct OutputConfiguration;

namespace output {

class Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	virtual void storeSettings(const Step &step) =0;
	virtual void storeSettings(size_t steps) =0;
	virtual void storeSettings(const std::vector<size_t> &steps) =0;

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution) =0;
	virtual void finalize() =0;

	virtual ~Store() {};

protected:
	Store(const OutputConfiguration &configuration): _configuration(configuration) {};

	const OutputConfiguration &_configuration;
};

}
}

#endif /* SRC_OUTPUT_STORE_H_ */
