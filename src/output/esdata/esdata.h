
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include "../store.h"
#include <cstdlib>

namespace espreso {
namespace output {

class Esdata: public MeshStore {

public:
	Esdata(const Mesh &mesh, const std::string &path): MeshStore(mesh, path) { };

	void store(double shrinkSubdomain, double shrinkCluster);

private:
	void coordinates(const Coordinates &coordinates);
	void elements(const Mesh &mesh);
	void materials(const Mesh &mesh, const std::vector<Material> &materials);
	void boundaryConditions(const Coordinates &coordinates, const std::vector<BoundaryCondition*> &conditions, size_t DOFs);
	void boundaries(const Mesh &mesh);
};


}
}


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
