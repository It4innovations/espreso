
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace espreso {
namespace input {

class API: public Loader {

public:
	static void load(
			APIMesh &mesh,
			const std::vector<eslocal> &eType,
			std::vector<std::vector<eslocal> > &eNodes,
			std::vector<eslocal> &neighbours,
			size_t size, const esglobal *ids)
	{
		ESINFO(OVERVIEW) << "Set mesh through API";
		API api(mesh, eType, eNodes, neighbours, size, ids);
		api.fill();
	}

protected:
	API(
			APIMesh &mesh,
			const std::vector<eslocal> &eType,
			std::vector<std::vector<eslocal> > &eNodes,
			std::vector<eslocal> &neighbours,
			size_t size, const eslocal *ids)
	: Loader(mesh), _mesh(mesh), _eType(eType), _eNodes(eNodes), _neighbours(neighbours), _size(size), _ids(ids) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials) { }; // unimportant
	void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes) { ESINFO(GLOBAL_ERROR) << "Implement settings for API."; }
	void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);

	void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) { }

private:
	APIMesh &_mesh;
	const std::vector<eslocal> &_eType;
	std::vector<std::vector<eslocal> > &_eNodes;
	std::vector<eslocal> &_neighbours;
	size_t _size;
	const eslocal *_ids;
};

}
}




#endif /* INPUT_API_API_H_ */
