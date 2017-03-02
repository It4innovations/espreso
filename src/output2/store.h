
#ifndef SRC_OUTPUT_STORE_H_
#define SRC_OUTPUT_STORE_H_

#include <string>
#include <vector>

namespace espreso {

class Mesh;
struct OutputConfiguration;
struct Step;

namespace output {

class Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	virtual void storeSettings(const Step &step) =0;
	virtual void finalize() {};

	virtual ~Store() {};

protected:
	Store(const OutputConfiguration &configuration, const Mesh *mesh, const std::string &path)
	:_configuration(configuration), _mesh(mesh), _path(path) {};

	const OutputConfiguration &_configuration;
	const Mesh *_mesh;
	std::string _path;

};

}
}



#endif /* SRC_OUTPUT_STORE_H_ */
