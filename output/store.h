
#ifndef OUTPUT_STORE_H_
#define OUTPUT_STORE_H_

#include "esmesh.h"

namespace esoutput {


class MeshStore {

public:
	virtual void store(double shrinkSubdomain, double shringCluster) = 0;

	virtual ~MeshStore() {};

protected:
	MeshStore(const mesh::Mesh &mesh, const std::string &path): _mesh(mesh), _path(path) {};

	const mesh::Mesh &_mesh;
	std::string _path;
};

class ResultStore: public MeshStore {

public:
	virtual void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster) = 0;

	virtual ~ResultStore() {};

protected:
	ResultStore(const mesh::Mesh &mesh, const std::string &path): MeshStore(mesh, path) {};
};

template <class TStore>
class Store {

public:
	Store(const mesh::Mesh &mesh, const std::string &path): _store(mesh, path) { };

	void store(double shrinkSubdomain = 1, double shringCluster = 1)
	{
		_store.store(shrinkSubdomain, shringCluster);
	}

	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain = 1, double shringCluster = 1)
	{
		_store.store(displacement, shrinkSubdomain, shringCluster);
	}

private:
	TStore _store;
};

}



#endif /* OUTPUT_STORE_H_ */
