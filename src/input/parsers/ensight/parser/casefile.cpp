
#include "casefile.h"
#include "basis/utilities/sysutils.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "basis/io/inputfile.h"

#include <cstring>
#include <cstdio>

using namespace espreso;

EnsightCasefile::EnsightCasefile(const std::string &path)
: path(path), type(Type::Ensight_Gold)
{
	parse();
}

void EnsightCasefile::parse()
{
	// currently read only geometry
	Metadata casefile;
	casefile.read(path);

	const char* current = casefile.begin;
	while (current < casefile.end) {
		while (current < casefile.end && *current++ != '\n'); // start at new line
		if (memcmp(current, "model:", 6) == 0) {
			while (*current++ != ':');
			while (*current == ' ' || *current == '\t') { ++current; }
			const char* begin = current;
			while (*current != '\n' && *current != '\r') { ++current; };
			geometry = std::string(begin, current);
			break;
		}
	}

	geometry = utils::getFileDirectory(path) + "/" + geometry;
}
