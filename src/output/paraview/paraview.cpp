
#include "paraview.h"

using namespace espreso::store;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Store(mesh, path)
{
	// constructor
	// save mesh to the memory
}

void Paraview::store(std::vector<std::vector<double> > &displasment, double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE RESULT TO PARAVIEW\n";
}
