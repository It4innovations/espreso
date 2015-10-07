
#ifndef OUTPUT_STORE_H_
#define OUTPUT_STORE_H_

#include "esmesh.h"

namespace esoutput {


class MeshStore {

public:
	virtual void store(const mesh::Mesh &mesh) = 0;

	virtual ~MeshStore() {};
};

template <class TStore>
class Store {

public:
	Store(const std::string &path, int rank, int size): _store(path, rank, size) { };

	void store(const mesh::Mesh &mesh)
	{
		_store.store(mesh);
	}

private:
	TStore _store;
};

}



#endif /* OUTPUT_STORE_H_ */
