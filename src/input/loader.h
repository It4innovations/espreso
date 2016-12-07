
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esbasis.h"
#include "esmesh.h"
#include <string>

namespace espreso {

struct GlobalConfiguration;

namespace input {

class Loader {

public:
	static void load(const GlobalConfiguration &configuration, Mesh &mesh, size_t index, size_t size);

	void fill();

	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges) = 0;
	virtual void materials(std::vector<Material> &materials) = 0;
	virtual void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours) = 0;
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes) = 0;

	virtual void open() {};
	virtual void close() {};

	virtual bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners) = 0;

	void boundaryConditions();

protected:
	Loader(const GlobalConfiguration &configuration, Mesh &mesh): configuration(configuration), mesh(mesh)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &environment->MPIrank);
		MPI_Comm_size(MPI_COMM_WORLD, &environment->MPIsize);
	}
	virtual ~Loader() {};

	const GlobalConfiguration &configuration;
	Mesh &mesh;
};

}
}



#endif /* INPUT_LOADER_H_ */
