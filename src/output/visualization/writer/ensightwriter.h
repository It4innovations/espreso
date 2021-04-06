
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_

#include "basis/io/outputfile.h"
#include "mesh/element.h"

namespace espreso {

struct EnsightWriter {
	static std::string codetotype(int code)
	{
		switch (static_cast<Element::CODE>(code)) {

		case Element::CODE::POINT1: return "point";

		case Element::CODE::LINE2: return "bar2";

		case Element::CODE::TRIANGLE3: return "tria3";
		case Element::CODE::SQUARE4: return "quad4";

		case Element::CODE::TETRA4: return "tetra4";
		case Element::CODE::PYRAMID5: return "pyramid5";
		case Element::CODE::PRISMA6: return "penta6";
		case Element::CODE::HEXA8: return "hexa8";

		case Element::CODE::LINE3: return "bar3";

		case Element::CODE::TRIANGLE6: return "tria6";
		case Element::CODE::SQUARE8: return "quad8";

		case Element::CODE::TETRA10: return "tetra10";
		case Element::CODE::PYRAMID13: return "pyramid13";
		case Element::CODE::PRISMA15: return "penta15";
		case Element::CODE::HEXA20: return "hexa20";

		default:
			return "";
		}
	}
};

struct EnsightBinaryWriter: public OutputFilePack {

	void format()
	{
		description("C Binary");
	}

	void description(const std::string &description)
	{
		insert(description.size(), description.data());
		insert(80 - description.size(), '\0');
	}

	void int32(int value)
	{
		insert(sizeof(int), &value);
	}

	void float32(float value)
	{
		insert(sizeof(float), &value);
	}

	void enode(int value)
	{
		insert(sizeof(int), &value);
	}

	void eend()
	{

	}
};

struct EnsightASCIIWriter: public OutputFilePack {

	void format()
	{

	}

	void description(const std::string &description)
	{
		insert(description.size(), description.data());
		insert(1, '\n');
	}

	void int32(int value)
	{
		insert(snprintf(buffer, bsize, "%10d\n", value));
	}

	void float32(float value)
	{
		insert(snprintf(buffer, bsize, "%12.5e\n", value));
	}

	void enode(int value)
	{
		insert(snprintf(buffer, bsize, "%10d", value));
	}

	void eend()
	{
		insert(1, '\n');
	}
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_ */
