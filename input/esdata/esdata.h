
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"

namespace esinput {

class Esdata: public ExternalLoader {

public:
	Esdata(const Options &options, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void faces(mesh::Faces &faces) {};
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries, std::vector<int> &neighbours);

	void open() {};
	void close() {};

private:
	std::string _path;
	int _rank;
	int _size;
};

}


#endif /* INPUT_ESDATA_ESDATA_H_ */
