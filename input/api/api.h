
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace esinput {

class API: public APILoader {

public:
	// TODO: elements with various DOFS
	API(std::vector<std::vector<eslocal> > &eIndices, std::vector<eslocal> &neighbours, size_t size, esglobal *ids)
	: _DOFs(3), _eIndices(eIndices), _neighbours(neighbours), _size(size), _ids(ids) { };

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries, std::vector<int> &neighbours);

private:
	size_t _DOFs;
	std::vector<std::vector<eslocal> > &_eIndices;
	std::vector<eslocal> &_neighbours;
	size_t _size;
	esglobal *_ids;
};

}




#endif /* INPUT_API_API_H_ */
