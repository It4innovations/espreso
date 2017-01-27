
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"


namespace espreso {

class APIMesh;
struct ESPRESOInput;

namespace input {

class API {

public:
	static void load(
			const ESPRESOInput &configuration,
			APIMesh &mesh,
			eslocal indexBase,
			const std::vector<eslocal> &eType,
			std::vector<std::vector<eslocal> > &eNodes,
			std::vector<std::vector<eslocal> > &eDOFs,
			std::vector<std::vector<double> > &eMatrices,
			eslocal dirichletSize,
			eslocal *dirichletIndices,
			double *dirichletValues,
			std::vector<int> &neighbours,
			size_t size, const eslocal *l2g);

protected:
	API(const ESPRESOInput &configuration, APIMesh &mesh, eslocal offset): _configuration(configuration), _mesh(mesh), _offset(offset) {};

	void points(const std::vector<std::vector<eslocal> > &eNodes, size_t DOFsSize);
	void elements(const std::vector<eslocal> &eType, std::vector<std::vector<eslocal> > &eNodes, std::vector<std::vector<eslocal> > &eDOFs, std::vector<std::vector<double> > &eMatrices);
	void dirichlet(size_t dirichletSize, eslocal *dirichletIndices, double *dirichletValues);
	void clusterBoundaries(std::vector<int> &neighbours, size_t size, const eslocal *l2g);

	const ESPRESOInput &_configuration;
	APIMesh &_mesh;
	eslocal _offset;
};

}
}




#endif /* INPUT_API_API_H_ */
