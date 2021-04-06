
#ifndef SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_
#define SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_

#include "parser.h"

namespace espreso {

struct BlockFinish: public AbaqusParser {
	static size_t nSize;
	static size_t ncSize;

	static const char* asterik;
	static const char* comment;

	BlockFinish();
	BlockFinish& parse(const char* begin);
};

}



#endif /* SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_ */


/*

#ifndef SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_
#define SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_

#include "parser.h"

namespace espreso {

struct BlockEnd: public AbaqusParser {
	static size_t nSize;
	static size_t unixSize;
	static size_t winSize;

	static const char* nUpper;
	static const char* nLower;
	static const char* unixEnd;
	static const char* winEnd;

	BlockEnd();
	BlockEnd& parse(const char* begin);
};

}



#endif  SRC_INPUT_ABAQUS_PARSER_BLOCKEND_H_ */
