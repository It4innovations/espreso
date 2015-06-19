#include "loader.h"

size_t Loader::getLinesCount(const char* fileName)
{
	std::ifstream file(fileName);
	if (file.is_open())
	{
		TestEOL test;
		size_t size = count_if(std::istreambuf_iterator<char>(file),
				std::istreambuf_iterator<char>(), test);

		file.close();
		return size;
	}
	else
	{
		fprintf(stderr, "Cannot load file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

std::ostream& operator<<(std::ostream& os, const Ansys &a)
{
	os << "Coordinates: " << a.coordinates() << "\n";
	os << "Elements: " << a.elements() << "\n";
	os << "Properties: \n";
	std::map<std::string, std::string>::const_iterator it =
			a.coordinatesProperties().begin();
	for (; it != a.coordinatesProperties().end(); ++it)
	{
		os << it->first << " - " << it->second << "\n";
	}
	return os;
}
