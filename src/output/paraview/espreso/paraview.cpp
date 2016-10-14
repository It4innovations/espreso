
// Dummy Paraview file
// Do nothing

#include "../paraview.h"

using namespace espreso::output;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Store(mesh, path)
{
	ESINFO(GLOBAL_ERROR) << "Re-compile ESPRESO with Paraview support.";
}

void Paraview::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster) { }

void Paraview::storeGeometry(size_t timeStep) { }
void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) { }
void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) { }
void Paraview::finalize() { }
