
#ifndef SRC_INPUT_ABAQUS_PARSER_SSECTION_H_
#define SRC_INPUT_ABAQUS_PARSER_SSECTION_H_

#include "parser.h"
#include <unordered_map>

namespace espreso {

struct Entry {
	std::string element;
	int number;

};

struct SSection: public AbaqusParser {
	static size_t size;
	static const char* upper;
	static const char* lower;
	static const char* sentence;

	char Elset[MAX_NAME_SIZE];
	char Material[MAX_NAME_SIZE];

	SSection();
	SSection& parse(const char* begin);

protected:
	std::unordered_map<std::string, std::string> elset_mat_dict;
	std::unordered_map<std::string, int> mat_index_dict;
	std::unordered_map<std::string, Entry> elset_data_dict;

};

}


#endif /* SRC_INPUT_ABAQUS_COMMAND_H_ */
