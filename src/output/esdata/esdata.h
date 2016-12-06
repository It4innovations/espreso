
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include <cstdlib>

#include "esmesh.h"

namespace espreso {
namespace store {

class Esdata {

public:
	static void mesh(const Mesh &mesh, const std::string &path);

protected:
	Esdata(const Mesh &mesh, const std::string &path);

	void coordinates(const Coordinates &coordinates);
	void elements(const Mesh &mesh);
	void materials(const Mesh &mesh, const std::vector<Material> &materials);
	void settings(const Mesh &mesh);
	void boundaries(const Mesh &mesh);

	const Mesh &_mesh;
	const std::string _path;

	// sorted faces and edges
};


}
}


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
