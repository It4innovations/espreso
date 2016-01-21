
#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

int main(int argc, char** argv)
{
	if (argc < 4) {
		std::cerr << "Specify parameters: ANSYS_FILE NUMBER_OF_PARTS OUTPUT_LOCATION\n";
		exit(EXIT_FAILURE);
	}

	esconfig::mesh::subdomains = atoi(argv[2]);
	esconfig::mesh::fixPoints = 0;
	mesh::Mesh m;
	esinput::AnsysWorkbench loader(argc, argv, 0, 1);
	loader.load(m);
	esoutput::Esdata data(m, argv[3]);
	data.store(1, 1);
}


