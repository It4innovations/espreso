
#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

int main(int argc, char** argv)
{
	if (argc < 4) {
		std::cerr << "Specify parameters: ANSYS_FILE NUMBER_OF_PARTS OUTPUT_LOCATION\n";
		exit(EXIT_FAILURE);
	}

	mesh::Mesh m(0, 0);
	m.load(mesh::ANSYS, argc, argv);
	m.partitiate(atoi(argv[2]), 0);
	m.store(mesh::ESPRESO_OUTPUT, argv[3]);
}


