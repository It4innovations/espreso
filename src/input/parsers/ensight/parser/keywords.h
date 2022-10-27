
#ifndef SRC_INPUT_PARSERS_ENSIGHT_PARSER_KEYWORDS_H_
#define SRC_INPUT_PARSERS_ENSIGHT_PARSER_KEYWORDS_H_

#include "esinfo/eslog.h"
#include "input/parsers/distributedscanner.h"
#include "mesh/mesh.h"
#include "mesh/element.h"

#include <cstddef>
#include <vector>

namespace espreso {

struct EnsightKeywords {
	enum class Format {
		BINARY,
		ASCII,
		UNKNOWN
	};
	enum class IDs {
		OFF,
		GIVEN,
		ASSIGN,
		INGNORE,
		UNKNOWN
	};

	struct Header {
		Format format;
		IDs nodeIDs, elementIDs;

		Header(): format(Format::UNKNOWN), nodeIDs(IDs::UNKNOWN), elementIDs(IDs::UNKNOWN) {}
	};

	struct Keyword {
		size_t fileindex, offset;

		Keyword(): fileindex(0), offset(0) {}
		Keyword(size_t fileindex, size_t offset): fileindex(fileindex), offset(offset) {}
	};

	struct Part: Keyword {
		Part(): number(0) {}
		Part(size_t fileindex, size_t offset, int number, const char *name): Keyword(fileindex, offset), number(number)
		{
			memcpy(this->name, name, 80);
		}

		std::string getName()
		{
			auto end = name;
			while (*end != 0 && *end != ' ' && *end != '\n') { ++end; }
			return std::string(name, end);
		}

		int number;
		char name[80];
	};

	struct Coordinates: Keyword {
		Coordinates(): nn(0) {}
		Coordinates(size_t fileindex, size_t offset, int nn): Keyword(fileindex, offset), nn(nn) {}

		int nn;
	};

	struct Elements: Keyword {
		enum class Type {
			POINT,
			BAR2, BAR3,
			TRIA3, TRIA6, QUAD4, QUAD8,
			TETRA4, TETRA10, PYRAMID5, PYRAMID13, PENTA6, PENTA15, HEXA8, HEXA20,
			NSIDED, NFACED
		};

		Elements(): type(Type::POINT), ne(0) {}
		Elements(size_t fileindex, size_t offset, Type type, int ne): Keyword(fileindex, offset), type(type), ne(ne) {}

		Element::CODE getCode() const
		{
			Element::CODE code;
			switch (type) {
			case Elements::Type::POINT    : code = Element::CODE::POINT1   ; break;
			case Elements::Type::BAR2     : code = Element::CODE::LINE2    ; break;
			case Elements::Type::BAR3     : code = Element::CODE::LINE3    ; break;
			case Elements::Type::TRIA3    : code = Element::CODE::TRIANGLE3; break;
			case Elements::Type::TRIA6    : code = Element::CODE::TRIANGLE6; break;
			case Elements::Type::QUAD4    : code = Element::CODE::SQUARE4  ; break;
			case Elements::Type::QUAD8    : code = Element::CODE::SQUARE8  ; break;
			case Elements::Type::TETRA4   : code = Element::CODE::TETRA4   ; break;
			case Elements::Type::TETRA10  : code = Element::CODE::TETRA10  ; break;
			case Elements::Type::PYRAMID5 : code = Element::CODE::PYRAMID5 ; break;
			case Elements::Type::PYRAMID13: code = Element::CODE::PYRAMID13; break;
			case Elements::Type::PENTA6   : code = Element::CODE::PRISMA6  ; break;
			case Elements::Type::PENTA15  : code = Element::CODE::PRISMA15 ; break;
			case Elements::Type::HEXA8    : code = Element::CODE::HEXA8    ; break;
			case Elements::Type::HEXA20   : code = Element::CODE::HEXA20   ; break;
			case Elements::Type::NSIDED   : code = Element::CODE::SIZE     ; break; // not supported
			case Elements::Type::NFACED   : code = Element::CODE::SIZE     ; break; // not supported
			default: code = Element::CODE::SIZE;
			}
			return code;
		}

		int getSize() const
		{
			return Mesh::edata[(int)getCode()].nodes;
		}

		Type type;
		int ne;
	};

	Header header;
	std::vector<Part> parts;
	std::vector<Coordinates> coordinates;
	std::vector<Elements> elements;
};

template <typename Parser>
void fillScanner(const FilePack &file, DistributedScanner &scanner, EnsightKeywords &keywords)
{
	scanner.add("part"       , [&] (const char *c) { Parser::addPart(file, c, keywords); });
	scanner.add("coordinates", [&] (const char *c) { Parser::addCoordinates(file, c, keywords); }, [&] (const char *c) { return Parser::skipCoordinates(c); });

	scanner.add("point"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::POINT); }    , [&] (const char *c) { return Parser::skipElements(c,  1); });
	scanner.add("bar2"       , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::BAR2); }     , [&] (const char *c) { return Parser::skipElements(c,  2); });
	scanner.add("bar3"       , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::BAR3); }     , [&] (const char *c) { return Parser::skipElements(c,  3); });
	scanner.add("tria3"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::TRIA3); }    , [&] (const char *c) { return Parser::skipElements(c,  3); });
	scanner.add("tria6"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::TRIA6); }    , [&] (const char *c) { return Parser::skipElements(c,  6); });
	scanner.add("quad4"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::QUAD4); }    , [&] (const char *c) { return Parser::skipElements(c,  4); });
	scanner.add("quad8"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::QUAD8); }    , [&] (const char *c) { return Parser::skipElements(c,  8); });
	scanner.add("tetra4"     , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::TETRA4); }   , [&] (const char *c) { return Parser::skipElements(c,  4); });
	scanner.add("tetra10"    , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::TETRA10); }  , [&] (const char *c) { return Parser::skipElements(c, 10); });
	scanner.add("pyramid5"   , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::PYRAMID5); } , [&] (const char *c) { return Parser::skipElements(c,  5); });
	scanner.add("pyramid13"  , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::PYRAMID13); }, [&] (const char *c) { return Parser::skipElements(c, 13); });
	scanner.add("penta6"     , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::PENTA6); }   , [&] (const char *c) { return Parser::skipElements(c,  6); });
	scanner.add("penta15"    , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::PENTA15); }  , [&] (const char *c) { return Parser::skipElements(c, 15); });
	scanner.add("hexa8"      , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::HEXA8); }    , [&] (const char *c) { return Parser::skipElements(c,  8); });
	scanner.add("hexa20"     , [&] (const char *c) { Parser::addElements(file, c, keywords, EnsightKeywords::Elements::Type::HEXA20); }   , [&] (const char *c) { return Parser::skipElements(c, 20); });

	scanner.add("nsided"     , [&] (const char *c) { eslog::error("Ensight scanner error: ESPRESO does not support nsided elements.\n"); });
	scanner.add("nfaced"     , [&] (const char *c) { eslog::error("Ensight scanner error: ESPRESO does not support nfaced elements.\n"); });
}

struct EnsightASCIIGeometryKeywordParser {
	static void addPart(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type);

	static size_t skipCoordinates(const char *c);
	static size_t skipElements(const char *c, int enodes);
};

struct EnsightBinaryGeometryKeywordParser {
	static void addPart(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type);

	static size_t skipCoordinates(const char *c);
	static size_t skipElements(const char *c, int enodes);
};

struct EnsightASCIIVariableKeywordParser {
	static void addPart(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type);

	static size_t skipCoordinates(const char *c);
	static size_t skipElements(const char *c, int enodes);
};

struct EnsightBinaryVariableKeywordParser {
	static void addPart(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addCoordinates(const FilePack &file, const char *c, EnsightKeywords &keywords);
	static void addElements(const FilePack &file, const char *c, EnsightKeywords &keywords, EnsightKeywords::Elements::Type type);

	static size_t skipCoordinates(const char *c);
	static size_t skipElements(const char *c, int enodes);
};

}

#endif /* SRC_INPUT_PARSERS_ENSIGHT_PARSER_KEYWORDS_H_ */
