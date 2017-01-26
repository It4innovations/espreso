
#include "paraview.h"

using namespace espreso::store;

Paraview::Paraview(const OutputConfiguration &output, const Mesh &mesh, const std::string &path): Store(output, mesh, path)
{
	// constructor
	// save mesh to the memory
}

void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{

}

void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{

}

void Paraview::store(std::vector<std::vector<double> > &displasment)
{
	std::cout << "SAVE RESULT TO PARAVIEW\n";
}
