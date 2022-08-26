
#include "polyBoundaryMesh.h"

#include <cstring>

using namespace espreso;

OpenFOAMPolyBoundaryMesh::OpenFOAMPolyBoundaryMesh(const std::string &file)
{
	std::ifstream is(file);
	FoamFileHeader header; header.read(is);

	auto get = [&is] (std::string &name, std::string &value) {
		is >> name;
		while (is.peek() == ' ') is.get();
		std::getline(is, value);
	};

	int bcount;
	is >> bcount;
	is.ignore(256, '(');
	boundaries.resize(bcount);
	for (int b = 0; b < bcount; ++b) {
		is >> boundaries[b].name;
		is.ignore(256, '{');
		std::string name, value;
		get(name, value);
		if (memcmp(name.c_str(), "type", 4) == 0) {
			if (memcmp(value.c_str(), "patch", 5) == 0) {
				boundaries[b].type = OpenFOAMBoundaryType::patch;
			}
			if (memcmp(value.c_str(), "wall", 4) == 0) {
				boundaries[b].type = OpenFOAMBoundaryType::wall;
			}
			if (memcmp(value.c_str(), "processor", 9) == 0) {
				boundaries[b].type = OpenFOAMBoundaryType::processor;
			}
		}
		if (boundaries[b].type == OpenFOAMBoundaryType::wall || boundaries[b].type == OpenFOAMBoundaryType::processor) {
			get(name, value); // inGroup
		}
		is >> name >> boundaries[b].nFaces; is.ignore(256, '\n');
		is >> name >> boundaries[b].startFace; is.ignore(256, '\n');

		if (boundaries[b].type == OpenFOAMBoundaryType::processor) {
			get(name, value); // match tolerance
			get(name, value); // transform
			get(name, value); // myProcNo
			is >> name >> boundaries[b].neighbor; is.ignore(256, '\n');
		}
		is.ignore(256, '}');
	}
}


