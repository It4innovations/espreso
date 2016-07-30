
#include "paraview.h"

using namespace espreso::output;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Results(mesh, path)
{
	// constructor
	// save mesh to the memory
}

void Paraview::store(std::vector<std::vector<double> > &displasment, double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE RESULT TO PARAVIEW\n";
}
