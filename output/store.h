
#ifndef OUTPUT_STORE_H_
#define OUTPUT_STORE_H_

#include "esmesh.h"

namespace esoutput {


class MeshStore {

public:
	virtual void store(double shrinkSubdomain, double shringCluster) = 0;

	virtual ~MeshStore() {};

protected:
	MeshStore(const mesh::Mesh &mesh, const std::string path): _mesh(mesh), _path(path) {};

	const mesh::Mesh &_mesh;
	std::string _path;
};

class ResultStore: public MeshStore {

public:
	virtual void store(std::vector<std::vector<double> > &displacement, size_t dofs, double shrinkSubdomain, double shringCluster) = 0;

	virtual ~ResultStore() {};

protected:
	ResultStore(const mesh::Mesh &mesh, const std::string &path): MeshStore(mesh, path) {};
};

}



#endif /* OUTPUT_STORE_H_ */
