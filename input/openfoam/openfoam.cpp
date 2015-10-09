
#include "openfoam.h"

using namespace esinput;

OpenFOAM::OpenFOAM(int argc, char** argv, int rank, int size)
{
	if (argc < 2) {
		if (rank == 0) {
			std::cerr << "Specify the path to an example as the first command line attribute.\n";
		}
		exit(EXIT_FAILURE);
	}
	_path = argv[1];
}

void OpenFOAM::points(mesh::Coordinates &coordinates)
{

}


void OpenFOAM::elements(std::vector<mesh::Element*> &elements)
{

}




