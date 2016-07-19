
// Dummy VTK file
// Always store VTK Legacy format

#include "vtk.h"
#include "../vtk/full/vtk.h"

using namespace espreso::output;

Generic::Generic(const Mesh &mesh, const std::string &path): ResultStore(mesh, path)
{
	ESINFO(ALWAYS) << "Warning: ESPRESO not contains a library for saving generic VTK format. VTK legacy format is used.";
}

void Generic::store(double shrinkSubdomain, double shrinkCluster)
{
	VTK_Full output(_mesh, _path);
	output.store(shrinkSubdomain, shrinkCluster);
}

void Generic::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster)
{
	VTK_Full output(_mesh, _path);
	output.store(displacement, shrinkSubdomain, shrinkCluster);

}
