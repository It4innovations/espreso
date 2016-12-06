
// Dummy Paraview file
// Do nothing

#include "paraview.h"

using namespace espreso::store;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Store(mesh, path)
{
	ESINFO(GLOBAL_ERROR) << "Re-compile ESPRESO with Paraview support.";
}

void Paraview::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster) { }
