
#ifndef SRC_OLD_OUTPUT_DATASTORE_H_
#define SRC_OLD_OUTPUT_DATASTORE_H_

#include <string>
#include <vector>

namespace espreso {

class Mesh;

namespace store {

class DataStore {

protected:
	DataStore(const Mesh &mesh, const std::string &path): _mesh(mesh), _path(path) {};

	const Mesh &_mesh;
	std::string _path;
};

}
}



#endif /* SRC_OLD_OUTPUT_DATASTORE_H_ */
