
// Dummy Paraview file
// Do nothing

#include "paraview.h"
#include "../../basis/logging/logging.h"

using namespace espreso::store;

Paraview::Paraview(const OutputConfiguration &output, const Mesh &mesh, const std::string &path): Store(output, mesh, path)
{
	ESINFO(GLOBAL_ERROR) << "Re-compile ESPRESO with Paraview support.";
}

void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	ESINFO(GLOBAL_ERROR) << "Implement store property";
}

void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	ESINFO(GLOBAL_ERROR) << "Implement store property";
}

void Paraview::store(std::vector<std::vector<double> > &displacement) { }
