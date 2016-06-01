
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"
#include "esbasis.h"

namespace espreso {
namespace input {

class Esdata: public Loader {

public:
	static void load(Mesh &mesh, const Configuration &configuration, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from ESPRESO binary format from directory " << configuration.path;
		Esdata esdata(mesh, configuration, rank, size);
		esdata.fill();
	}

protected:
	Esdata(Mesh &mesh, const Configuration &configuration, int rank, int size)
	: Loader(mesh), _path(configuration.path), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates, size_t &DOFs);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials);
	void boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

private:
	std::string _path;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESDATA_ESDATA_H_ */
