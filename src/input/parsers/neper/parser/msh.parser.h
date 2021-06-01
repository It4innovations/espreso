
#ifndef SRC_INPUT_PARSERS_NEPER_PARSER_MSH_PARSER_H_
#define SRC_INPUT_PARSERS_NEPER_PARSER_MSH_PARSER_H_

#include "basis/io/inputfile.h"

#define NEPER_MAX_NAME_SIZE 50

namespace espreso {

struct MeshBuilder;

class NeperMshMesh {
	struct Keyword {
		size_t offset, begin, end;
		int rbegin, rend;

		Keyword(): offset(-1), begin(-1), end(-1), rbegin(-1), rend(-1) {}
		Keyword(InputFilePack &pack, const char *c);
	};

	struct KeywordWithSize: public Keyword {
		size_t size;

		KeywordWithSize(): size(0) {}
		KeywordWithSize(InputFilePack &pack, const char *c, const std::string &open, const std::string &close);
	};

	struct Format: public Keyword {
		static constexpr const char *open = "$MeshFormat", *close = "$EndMeshFormat";

		enum File: int {
			ASCII = 0,
			binary = 1
		};

		double version;
		int file_type, data_size;

		Format(): version(-1), file_type(-1), data_size(-1) {}
		Format(InputFilePack &pack, const char *c);
	};

	struct Nodes: public KeywordWithSize {
		static constexpr const char *open = "$Nodes", *close = "$EndNodes";

		Nodes() {}
		Nodes(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct Elements: public KeywordWithSize {
		static constexpr const char *open = "$Elements", *close = "$EndElements";

		Elements() {}
		Elements(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct Periodicity: public KeywordWithSize {
		static constexpr const char *open = "$Periodicity", *close = "$EndPeriodicity";

		Periodicity() {}
		Periodicity(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct NSets: public KeywordWithSize {
		static constexpr const char *open = "$NSets", *close = "$EndNSets";

		NSets() {}
		NSets(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct Fasets: public KeywordWithSize {
		static constexpr const char *open = "$Fasets", *close = "$EndFasets";

		Fasets() {}
		Fasets(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct NodePartitions: public KeywordWithSize {
		static constexpr const char *open = "$NodePartitions", *close = "$EndNodePartitions";

		NodePartitions() {}
		NodePartitions(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct PhysicalNames: public KeywordWithSize {
		static constexpr const char *open = "$PhysicalNames", *close = "$EndPhysicalNames";

		PhysicalNames() {}
		PhysicalNames(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct ElsetOrientations: public KeywordWithSize {
		static constexpr const char *open = "$ElsetOrientations", *close = "$EndElsetOrientations";

		ElsetOrientations() {}
		ElsetOrientations(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct ElementOrientations: public KeywordWithSize {
		static constexpr const char *open = "$ElementOrientations", *close = "$EndElementOrientations";

		ElementOrientations() {}
		ElementOrientations(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

	struct Groups: public KeywordWithSize {
		static constexpr const char *open = "$Groups", *close = "$EndGroups";

		Groups() {}
		Groups(InputFilePack &pack, const char *c): KeywordWithSize(pack, c, open, close) {}
	};

public:
	NeperMshMesh(InputFilePack &meshfile);

	void parse(MeshBuilder &mesh);

protected:
	void scan();
	void parseASCII(MeshBuilder &mesh);

	InputFilePack &_meshfile;
	std::vector<Format> _format;
	std::vector<Nodes> _nodes;
	std::vector<Elements> _elements;
	std::vector<Periodicity> _periodicity;
	std::vector<NSets> _nsets;
	std::vector<Fasets> _fasets;
	std::vector<NodePartitions> _nodePartitions;
	std::vector<PhysicalNames> _physicalNames;
	std::vector<ElsetOrientations> _elsetOrientations;
	std::vector<ElementOrientations> _elementOrientations;
	std::vector<Groups> _groups;
};
}





#endif /* SRC_INPUT_PARSERS_NEPER_PARSER_MSH_PARSER_H_ */
