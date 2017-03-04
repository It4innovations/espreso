
#ifndef SRC_OLD_OUTPUT_RESULTSTORE_H_
#define SRC_OLD_OUTPUT_RESULTSTORE_H_

#include <string>
#include <vector>

namespace espreso {

class Mesh;
class OutputConfiguration;
enum class Property;

namespace store {

enum class ElementType {
	NODES,
	EDGES,
	FACES,
	ELEMENTS
};

class ResultStore {

public:
	const OutputConfiguration& configuration() const { return _output; };

	virtual void storeGeometry(size_t timeStep = -1) = 0;
	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) = 0;
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) = 0;
	virtual void finalize() =0;

	virtual ~ResultStore() {};

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
	:_output(output), _mesh(mesh), _path(path) {};

	const OutputConfiguration &_output;
	const Mesh *_mesh;
	std::string _path;
};

}
}



#endif /* SRC_OLD_OUTPUT_RESULTSTORE_H_ */
