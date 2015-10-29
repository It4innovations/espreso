
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include "../store.h"
#include <cstdlib>

namespace esoutput {

class Esdata: public MeshStore {

public:
	Esdata(const mesh::Mesh &mesh, const std::string &path): MeshStore(mesh, path) { };

	void store(double shrinkSubdomain, double shrinkCluster);

private:
	void coordinates(const mesh::Coordinates &coordinates);
	void elements(const mesh::Mesh &mesh);
	void boundaryConditions(const mesh::Coordinates &coordinates);
	void boundaries(const mesh::Mesh &mesh);
};


};


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
