
#ifndef SRC_INPUT_ABAQUS_PARSER_EBLOCK_H_
#define SRC_INPUT_ABAQUS_PARSER_EBLOCK_H_

#include "parser.h"
#include <string>

namespace espreso {

struct MeshBuilder;
struct ET;

enum class ETYPE {
		D1LINE_2NODES,
		D1LINE_3NODES,

		D2SOLID_4NODES,
		D2SOLID_6NODES,
		D2SOLID_8NODES,

		D3SOLID_4NODES,
		D3SOLID_8NODES,
		D3SOLID_10NODES,
		D3SOLID_20NODES,

		D2SURFACE,
		D3SURFACE,

		UNIVERSAL,

		UNKNOWN
	};




struct EList: public AbaqusParser {
	static size_t size;
	static const char* upper;
	static const char* lower;
	static const char* sentence;
	//std::vector<std::string> TYPE;
	char TYPE[MAX_NAME_SIZE];
	esint NUM_NODES, Solkey, NDMAX, NDSEL;

	esint lineSize, elementSize, lineEndSize;
	esint indexSize, indexLength, valueSize, valueLength;

	esint NUM_NODES_LINE1,NUM_NODE_LINE2;
	EList();
	EList& parse(const char* begin);
	bool readData(MeshBuilder &mesh);

	void fixOffsets(std::vector<size_t> &dataOffsets);

protected:
	bool solid( MeshBuilder &mesh);
	/*	bool solidBlock();

	void fixOffsets(std::vector<size_t> &dataOffsets);


	bool boundary(const std::vector<ET> &et, PlainAbaqusData &mesh);*/
};

}


#endif /* SRC_INPUT_ABAQUS_PARSER_EBLOCK_H_ */
