
#include "keywords.h"

using namespace espreso;

void EnsightASCIIGeometryKeywordParser::addPart(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	while (*c++ != '\n');
	int number = std::atof(c);
	while (*c++ != '\n');
	keywords.parts.push_back(EnsightKeywords::Part(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), number, c));
}

void EnsightBinaryGeometryKeywordParser::addPart(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	int number;
	memcpy(&number, c + 80, sizeof(int));
	keywords.parts.push_back(EnsightKeywords::Part(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), number, c + 84));
}

void EnsightASCIIVariableKeywordParser::addPart(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	const char *_c = c;
	while (*c++ != '\n');
	int number = std::atof(c);
	while (*c++ != '\n');
	keywords.parts.push_back(EnsightKeywords::Part(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), number, _c));
}

void EnsightBinaryVariableKeywordParser::addPart(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	int number;
	memcpy(&number, c + 80, sizeof(int));
	keywords.parts.push_back(EnsightKeywords::Part(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), number, c));
}


void EnsightASCIIGeometryKeywordParser::addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	while (*c++ != '\n');
	int count = atof(c);
	while (*c++ != '\n');
	keywords.coordinates.push_back(EnsightKeywords::Coordinates(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), count));
}

size_t EnsightASCIIGeometryKeywordParser::skipCoordinates(const char *c)
{
	while (*c++ != '\n');
	int count = atof(c);
	return 13 * 3 * count;
}

void EnsightBinaryGeometryKeywordParser::addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	int count;
	memcpy(&count, c + 80, sizeof(int));
	keywords.coordinates.push_back(EnsightKeywords::Coordinates(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin) + 80 + sizeof(int), count));
}

size_t EnsightBinaryGeometryKeywordParser::skipCoordinates(const char *c)
{
	int count;
	memcpy(&count, c + 80, sizeof(int));
	return count * 3 * sizeof(float);
}

void EnsightASCIIVariableKeywordParser::addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	// we need to get count from geometry file
	while (*c++ != '\n');
	keywords.coordinates.push_back(EnsightKeywords::Coordinates(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), -1));
}

size_t EnsightASCIIVariableKeywordParser::skipCoordinates(const char *c)
{
	return 0;
}

void EnsightBinaryVariableKeywordParser::addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords)
{
	keywords.coordinates.push_back(EnsightKeywords::Coordinates(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin) + 80, -1));
}

size_t EnsightBinaryVariableKeywordParser::skipCoordinates(const char *c)
{
	return 0;
}


void EnsightASCIIGeometryKeywordParser::addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type)
{
	while (*c++ != '\n');
	int nn = atof(c);
	while (*c++ != '\n');
	keywords.elements.push_back(EnsightKeywords::Elements(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), type, nn));
}

size_t EnsightASCIIGeometryKeywordParser::skipElements(const char *c, int enodes)
{
	while (*c++ != '\n');
	int nn = atof(c);
	return (enodes * 10 + 1) * nn;
}

void EnsightBinaryGeometryKeywordParser::addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type)
{
	int nn;
	memcpy(&nn, c + 80, sizeof(int));
	keywords.elements.push_back(EnsightKeywords::Elements(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin) + 80 + sizeof(int), type, nn));
}

size_t EnsightBinaryGeometryKeywordParser::skipElements(const char *c, int enodes)
{
	int nn; memcpy(&nn, c + 80, sizeof(int));
	return enodes * nn * sizeof(int);
}

void EnsightASCIIVariableKeywordParser::addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type)
{
	while (*c++ != '\n');
	keywords.elements.push_back(EnsightKeywords::Elements(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin), type, -1));
}

size_t EnsightASCIIVariableKeywordParser::skipElements(const char *c, int enodes)
{
	return 0;
}

void EnsightBinaryVariableKeywordParser::addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type)
{
	keywords.elements.push_back(EnsightKeywords::Elements(file.fileindex, file.distribution[info::mpi::rank] + (c - file.begin) + 80, type, -1));
}

size_t EnsightBinaryVariableKeywordParser::skipElements(const char *c, int enodes)
{
	return 0;
}

