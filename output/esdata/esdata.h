
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include "../store.h"
#include <cstdlib>

namespace esoutput {

class Esdata: public MeshStore {

public:
	Esdata(const std::string &path, int rank, int size): _path(path), _rank(rank), _size(size) { };

	void store(const mesh::Mesh &mesh, double shrinkSubdomain, double shrinkCluster);

private:
	std::string _path;
	int _rank;
	int _size;
};


};


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
