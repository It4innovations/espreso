
#ifndef INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_
#define INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_

#include "../loader.h"

namespace espreso {

struct ESPRESOInput;

namespace input {

class ESPRESOBinaryFormat: public Loader {

public:
	static void load(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size);

protected:
	ESPRESOBinaryFormat(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size)
	: Loader(mesh), _esdata(configuration), _rank(rank), _size(size) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void materials(std::vector<Material*> &materials);
	void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges);
	bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

private:
	const ESPRESOInput &_esdata;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESPRESO_ESPRESOBINARYFORMAT_H_ */
