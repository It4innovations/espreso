
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_TOKENIZER_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_TOKENIZER_H_

#include <string>
#include <cstring>
#include <fstream>

namespace espreso {
namespace oftoken {

struct dictionary { struct { const char *begin, *end; } keyword, value; };

bool isEmpty(const char &c)
{
	return c == ' ' || c == '\n' || c == '\r' || c == '\t';
}

const char* toNonEmpty(const char* c)
{
	while (isEmpty(*c)) { ++c; }
	return c;
}

const char* toEmpty(const char* c)
{
	while (!isEmpty(*c)) { ++c; }
	return c;
}

const char* toEnd(const char* c)
{
	while (*c++ != ';');
	return c;
}

const char* toFoamFile(const char* c)
{
	while (true) {
		while (*c != 'F') { ++c; }
		if (memcmp(c, "FoamFile", 8) == 0) {
			while (*c++ != '{');
			return c;
		}
		++c;
	}
}

void toFoamFile(std::ifstream &is)
{
	while (is.good()) {
		char line[256];
		is.getline(line, 256);
		if (memcmp(line, "FoamFile", 8) == 0) {
			is.getline(line, 256);
			return;
		}
	}
}


dictionary dict(const char* c)
{
	dictionary dict;
	dict.keyword.begin = toNonEmpty(c);
	dict.keyword.end = toEmpty(dict.keyword.begin);
	dict.value.begin = toNonEmpty(dict.keyword.end);
	dict.value.end = toEnd(dict.value.begin);
	return dict;
}

}
}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_TOKENIZER_H_ */
