
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"
#include "esbasis.h"

#include "../../config/input.h"

namespace espreso {
namespace input {

class Esdata: public Loader {

public:
	static void load(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from ESPRESO binary format from directory " << configuration.path;
		Esdata esdata(configuration, mesh, rank, size);
		esdata.fill();
	}

protected:
	Esdata(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size)
	: Loader(mesh), _esdata(configuration), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges);
	bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
	{
		mesh.partitiate(_esdata.domains);
		return true;
	}

private:
	const ESPRESOInput &_esdata;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESDATA_ESDATA_H_ */
