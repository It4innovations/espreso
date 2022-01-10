
#include "casefile.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/sysutils.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "basis/io/inputfile.h"

#include <cstring>
#include <cstdio>
#include <string>

using namespace espreso;

EnsightCasefile::EnsightCasefile(const std::string &path)
: path(path), type(Type::Ensight_Gold)
{
	parse();
}

void EnsightCasefile::parse()
{
	// currently read only geometry
	Metadata casefile(path);

	bool ingeometry = false, invariable = false, intime = false;
	size_t timesteps = 0, timeset = 0;

	const char* current = casefile.begin;
	while (current < casefile.end) {
		while (current < casefile.end && *current++ != '\n'); // go to the line end
		while (current < casefile.end && *current == '\n') { current++; } // start at new line
		if (current + 8 < casefile.end && memcmp(current, "GEOMETRY", 8) == 0) {
			ingeometry = true;
			invariable = intime = false;
			continue;
		}
		if (current + 8 < casefile.end && memcmp(current, "VARIABLE", 8) == 0) {
			invariable = true;
			ingeometry = intime = false;
			continue;
		}
		if (current + 4 < casefile.end && memcmp(current, "TIME", 4) == 0) {
			intime = true;
			ingeometry = invariable = false;
			continue;
		}
		if (ingeometry && current + 6 < casefile.end && memcmp(current, "model:", 6) == 0) {
			while (*current++ != ':');
			while (*current == ' ' || *current == '\t') { ++current; }
			const char* begin = current;
			while (*current != '\n') { ++current; };
			geometry = std::string(begin, current);
		}
		if (invariable) {
			Variable::Type type;
			int dimension = 0, time = -1;
			if (current + 6 < casefile.end && memcmp(current, "scalar", 6) == 0) {
				dimension = 1;
			}
			if (current + 6 < casefile.end && memcmp(current, "vector", 6) == 0) {
				dimension = 3;
			}
			if (current + 15 < casefile.end && memcmp(current + 7, "per node", 8) == 0) {
				type = Variable::Type::NODE;
			}
			if (current + 18 < casefile.end && memcmp(current + 7, "per element", 11) == 0) {
				type = Variable::Type::ELEMENT;
			}
			if (dimension) {
				const char *start = current;
				while (*start++ != ':');
				while (*++start == ' ');
				const char *end = start;
				while (*++end != '\n');
				std::vector<std::string> keys = Parser::split(std::string(start, end), " ");
				switch (keys.size()) {
				case 3:
					time = std::stod(keys[0]) - 1;
					variables.push_back(Variable{dimension, type, time, keys[1], keys[2]});
					break;
				case 2:
					variables.push_back(Variable{dimension, type, time, keys[0], keys[1]});
					break;
				default:
					eslog::error("Unknown variable in the case file.\n");
				}
			}
		}
		if (intime) {
			if (current + 15 < casefile.end && memcmp(current, "number of steps", 15) == 0) {
				while (*current++ != ':');
				timesteps = strtol(current, NULL, 10);
			}
			if (current + 8 < casefile.end && memcmp(current, "time set", 8) == 0) {
				while (*current++ != ':');
				timeset = strtol(current, NULL, 10);
				if (times.size() < timeset) {
					times.resize(timeset);
				}
				--timeset;
			}
			if (current + 11 < casefile.end && memcmp(current, "time values", 11) == 0) {
				const char *start = current;
				while (*start++ != ':');
				for (size_t s = 0; s < timesteps; ++s) {
					char *end;
					times[timeset].push_back(strtod(start, &end));
					start = end;
				}
			}
		}
	}

	geometry = utils::getFileDirectory(path) + "/" + geometry;
}
