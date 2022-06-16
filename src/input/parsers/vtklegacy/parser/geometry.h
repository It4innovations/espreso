
#ifndef SRC_INPUT_FORMATS_VTKLEGACY_PARSER_GEOMETRY_H_
#define SRC_INPUT_FORMATS_VTKLEGACY_PARSER_GEOMETRY_H_

#include "basis/io/inputfile.h"

#include <string>
#include <vector>

namespace espreso {

struct MeshBuilder;

class VTKLegacyGeometry {
	enum class Format {
		BINARY,
		ASCII,
		UNKNOWN
	};

	enum class DataSet {
		UNSTRUCTURED_GRID,
		UNKNOWN
	};

	enum class DataSource {
		POINTS,
		CELLS
	};

	// empty constructors are only for unpacking
	struct Keyword {
		size_t fileindex;
		size_t offset, begin, end;
		int rank;

		Keyword(): fileindex((size_t)-1), offset((size_t)-1), begin((size_t)-1), end((size_t)-1), rank(-1) {}
		Keyword(InputFilePack &pack, const char *c);
	};

	struct Header: public Keyword {
		Format format;
		DataSet dataset;

		Header(): format(Format::UNKNOWN), dataset(DataSet::UNKNOWN) {}
		Header(InputFilePack &pack, const char *c);
	};

	struct Points: public Keyword {
		size_t nn;

		Points(): nn(0) {}
		Points(InputFilePack &pack, const char *c);
	};

	struct Cells: public Keyword {
		size_t ne, size;

		Cells(): ne(0), size(0) {}
		Cells(InputFilePack &pack, const char *c);
	};

	struct CellTypes: public Keyword {
		size_t ne;

		CellTypes(): ne(0) {}
		CellTypes(InputFilePack &pack, const char *c);
	};

	struct Data: public Keyword {
		DataSource source;

		Data(): source(DataSource::POINTS) {}
		Data(InputFilePack &pack, DataSource source, const char *c);
	};

public:
	VTKLegacyGeometry(InputFilePack &pack);

	void scan();
	void parse(MeshBuilder &mesh, const std::vector<std::string> &names);

protected:
	void header();
	void scanBinary();
	void scanASCII();
	void parseBinary(MeshBuilder &mesh, const std::vector<std::string> &names);
	void parseASCII(MeshBuilder &mesh, const std::vector<std::string> &names);

	InputFilePack &_pack;

	std::vector<Header> _header;
	std::vector<Points> _points;
	std::vector<Cells> _cells;
	std::vector<CellTypes> _cellTypes;
	std::vector<Data> _pointData, _cellData;
};
}



#endif /* SRC_INPUT_FORMATS_VTKLEGACY_PARSER_GEOMETRY_H_ */
