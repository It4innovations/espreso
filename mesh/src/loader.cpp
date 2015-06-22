
#include "loader.h"

using namespace mesh;

size_t Loader::getLinesCount(const char* fileName)
{
	std::ifstream file(fileName);
	if (file.is_open()) {
		TestEOL test;
		size_t size = count_if(
				std::istreambuf_iterator<char>(file),
				std::istreambuf_iterator<char>(),
				test);

		file.close();
		return size;
	} else {
		fprintf(stderr, "Cannot load file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Ansys &a)
{
	os << "Coordinates: " << a.coordinates() << "\n";
	os << "Elements: " << a.elements() << "\n";
	os << "Properties: \n";
	for (size_t i = 0; i < a._coordinatesProperty.size(); i++) {
		os << a._coordinatesProperty[i] << "\n";
	}
	return os;
}
