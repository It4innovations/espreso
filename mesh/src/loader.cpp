#include "loader.h"

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

