
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace esinput {

class API: public ExternalLoader {

public:
	API(std::vector<std::vector<eslocal> > eIndices, std::vector<std::vector<double> > eMatrix)
	: eIndices(eIndices), eMatrix(eMatrix) { };

	void points(mesh::Coordinates &coordinates) { /* empty in case of API */ };
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates) { /* empty in case of API */ };
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries) { /* empty in case of API */ };

	void open() {};
	void close() {};

private:
	std::vector<std::vector<eslocal> > eIndices;
	std::vector<std::vector<double> > eMatrix;
};

}




#endif /* INPUT_API_API_H_ */
