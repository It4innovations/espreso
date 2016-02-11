
#include "openfoam.h"

using namespace esinput;

OpenFOAM::OpenFOAM(const Options &options, int rank, int size)
{
	_path = options.path;
}

void OpenFOAM::points(mesh::Coordinates &coordinates)
{

}


void OpenFOAM::elements(std::vector<mesh::Element*> &elements)
{

}

void OpenFOAM::boundaryConditions(mesh::Coordinates &coordinates)
{

}

void OpenFOAM::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{

}




