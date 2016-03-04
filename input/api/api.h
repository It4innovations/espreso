
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace esinput {

class API: public APILoader {

public:
	// TODO: elements with various DOFS
	API(std::vector<std::vector<eslocal> > &eIndices, std::vector<eslocal> &neighbours, size_t size, esglobal *ids)
	: DOFs(3), eIndices(eIndices), neighbours(neighbours), size(size), ids(ids) { };

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

private:
	size_t DOFs;
	std::vector<std::vector<eslocal> > &eIndices;
	std::vector<eslocal> &neighbours;
	size_t size;
	esglobal *ids;
};

}




#endif /* INPUT_API_API_H_ */
