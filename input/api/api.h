
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace esinput {

class API: public APILoader {

public:
	// TODO: elements with various DOFS
	API(std::vector<std::vector<eslocal> > &eIndices): DOFS(3), eIndices(eIndices) { };

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);

private:
	size_t DOFS;
	std::vector<std::vector<eslocal> > &eIndices;
};

}




#endif /* INPUT_API_API_H_ */
