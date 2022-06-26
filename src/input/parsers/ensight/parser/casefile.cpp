
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
: path(path), type(Type::UNKNOWN)
{
	parse();
}

void EnsightCasefile::parse()
{
	// currently read only geometry
	Metadata casefile(path);

	bool informat = false, ingeometry = false, invariable = false, intime = false;
	size_t timesteps = 0, timeset = 0;

	const char* current = casefile.begin;
	auto check = [&current,&casefile] (const char *parameter) {
		return current + strlen(parameter) < casefile.end && memcmp(current, parameter, strlen(parameter)) == 0;
	};

	while (current < casefile.end) {
		while (current < casefile.end && *current++ != '\n'); // go to the line end
		while (current < casefile.end && *current == '\n') { current++; } // start at new line
		if (check("FORMAT")) {
			informat = true;
			ingeometry = invariable = intime = false;
			continue;
		}
		if (check("GEOMETRY")) {
			ingeometry = true;
			informat = invariable = intime = false;
			continue;
		}
		if (check("VARIABLE")) {
			invariable = true;
			informat = ingeometry = intime = false;
			continue;
		}
		if (check("TIME")) {
			intime = true;
			informat = ingeometry = invariable = false;
			continue;
		}
		if (informat && check("type:")) {
			while (*current++ != ':');
			while (*current == ' ' || *current == '\t') { ++current; }
			if (check("ensight gold")) {
				type = Type::ENSIGHT_GOLD;
			}
		}
		if (ingeometry && check("model:")) {
			while (*current++ != ':');
			while (*current == ' ' || *current == '\t') { ++current; }
			const char* begin = current;
			while (*current != '\n') { ++current; };
			geometry = std::string(begin, current);
		}
		if (invariable) {
			Variable::Type type;
			int dimension = 0, time = -1;
			if (check("scalar")) {
				dimension = 1;
			}
			if (check("vector")) {
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
			if (check("number of steps")) {
				while (*current++ != ':');
				timesteps = strtol(current, NULL, 10);
			}
			if (check("time set")) {
				while (*current++ != ':');
				timeset = strtol(current, NULL, 10);
				if (times.size() < timeset) {
					times.resize(timeset);
				}
				--timeset;
			}
			if (check("time values")) {
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

	if (type == Type::UNKNOWN) {
		eslog::globalerror("Only Ensight Gold is supported.\n");
	}

	geometry = utils::getFileDirectory(path) + "/" + geometry;
}
