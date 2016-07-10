
#include "esoutput.h"
#include "esmesh.h"

using namespace espreso::output;

Generic::Generic(const Mesh &mesh, const std::string &path): ResultStore(mesh, path)
{
	// constructor
}

void Generic::store(double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE GENERIC VTK DATA\n";
}

void Generic::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE GENERIC VTK RESULT\n";
}



