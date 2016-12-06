
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"
#include "esbasis.h"

#include "../../config/description.h"

namespace espreso {
namespace input {

class Esdata: public Loader {

public:
	static void load(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from ESPRESO binary format from directory " << configuration.esdata.path;
		Esdata esdata(configuration, mesh, rank, size);
		esdata.fill();
	}

protected:
	Esdata(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size)
	: Loader(configuration, mesh), _esdata(configuration.esdata), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void materials(std::vector<Material> &materials);
	void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours);

private:
	ESPRESOInput _esdata;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESDATA_ESDATA_H_ */
