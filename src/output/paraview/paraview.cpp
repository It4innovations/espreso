
#include "paraview.h"

using namespace espreso::store;

Paraview::Paraview(const OutputConfiguration &output, const Mesh &mesh, const std::string &path): Store(output, mesh, path)
{
	// constructor
	// save mesh to the memory
}

void Paraview::store(std::vector<std::vector<double> > &displasment)
{
	std::cout << "SAVE RESULT TO PARAVIEW\n";
}
