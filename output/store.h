
#ifndef OUTPUT_STORE_H_
#define OUTPUT_STORE_H_

#include "esmesh.h"

namespace esoutput {


class MeshStore {

public:
	virtual void store(const mesh::Mesh &mesh, double shrinkSubdomain, double shringCluster) = 0;

	virtual ~MeshStore() {};
};

class ResultStore: public MeshStore {

public:
	virtual void store(const mesh::Mesh &mesh, std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster) = 0;

	virtual ~ResultStore() {};
};

template <class TStore>
class Store {

public:
	Store(const std::string &path, int rank, int size): _store(path, rank, size) { };

	void store(const mesh::Mesh &mesh, double shrinkSubdomain = 1, double shringCluster = 1)
	{
		_store.store(mesh, shrinkSubdomain, shringCluster);
	}

	void store(const mesh::Mesh &mesh, std::vector<std::vector<double> > &displacement, double shrinkSubdomain = 1, double shringCluster = 1)
	{
		_store.store(mesh, displacement, shrinkSubdomain, shringCluster);
	}

private:
	TStore _store;
};

}



#endif /* OUTPUT_STORE_H_ */
