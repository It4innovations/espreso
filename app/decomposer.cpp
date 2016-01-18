
#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

int main(int argc, char** argv)
{
	if (argc < 4) {
		std::cerr << "Specify parameters: ANSYS_FILE NUMBER_OF_PARTS OUTPUT_LOCATION\n";
		exit(EXIT_FAILURE);
	}

	mesh::Mesh m;
	esinput::AnsysWorkbench loader(argc, argv, 0, 1);
	loader.load(m);
	m.partitiate(atoi(argv[2]));
	esoutput::Esdata data(m, argv[3]);
	data.store(1, 1);
}


