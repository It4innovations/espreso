
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"
#include "esbasis.h"

namespace espreso {
namespace input {

class Esdata: public Loader {

public:
	static void load(Mesh &mesh, const Options &options, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from ESPRESO binary format from directory " << options.path;
		Esdata esdata(mesh, options, rank, size);
		esdata.fill();
	}

protected:
	Esdata(Mesh &mesh, const Options &options, int rank, int size)
	: Loader(mesh), _path(options.path), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates, size_t &DOFs);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials);
	void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);

private:
	std::string _path;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESDATA_ESDATA_H_ */
