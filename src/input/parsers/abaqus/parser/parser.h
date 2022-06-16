
#ifndef SRC_INPUT_ABAQUS_PARSER_PARSER_H_
#define SRC_INPUT_ABAQUS_PARSER_PARSER_H_

#include <cstddef>
#include <vector>
#include <string>

#define MAX_NAME_SIZE 50 // upper bound on name size
#define MAX_COMMAND_SIZE 150 // upper bound on command size
#define MAX_LINE_SIZE 512 // upper bound on line size
#define MAX_LINE_STEP 4   // sometimes we need to red more lines to get full information

namespace espreso {

struct BlockFinish;

struct AbaqusParser {
	static size_t offset;
	static const char* begin;
	static const char* end;

	size_t header;
	size_t first, last;
	int fRank, lRank;

	AbaqusParser()
	: header((size_t)-1),
	  first((size_t)-1), last((size_t)-1), fRank(-1), lRank(-1) {}

	void fillIndices(const char* header, const char* data);
	void fillIndices(const char* header, const char* first, const char* last);

	void fillDistribution(std::vector<BlockFinish> &blocksFinishs, std::vector<size_t> &distribution);
	const char* getFirst() const;
	const char* getLast() const;

	std::string command() const;

	void print(const char* data);
};
}



#endif /* SRC_INPUT_ABAQUS_PARSER_PARSER_H_ */
