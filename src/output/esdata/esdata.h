
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include <cstdlib>

#include "esmesh.h"

namespace espreso {
namespace output {

class Esdata {

public:
	static void mesh(const Mesh &mesh, const std::string &path);

protected:
	Esdata(const Mesh &mesh, const std::string &path);

	void coordinates(const Coordinates &coordinates);
	void elements(const Mesh &mesh);
	void materials(const Mesh &mesh, const std::vector<Material> &materials);
	void boundaries(const Mesh &mesh);

	const Mesh &_mesh;
	const std::string _path;
};


}
}


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
