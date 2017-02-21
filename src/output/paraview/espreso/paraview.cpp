
// Dummy Paraview file
// Do nothing

#include "../paraview.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::store;

Paraview::Paraview(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	ESINFO(GLOBAL_ERROR) << "Re-compile ESPRESO with Paraview support.";
}

void Paraview::storeGeometry(size_t timeStep) { }
void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) { }
void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) { }
void Paraview::finalize() { }
