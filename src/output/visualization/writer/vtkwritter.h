
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_

#include "basis/io/outputfile.h"
#include "mesh/element.h"

namespace espreso {

struct VTKASCIIWritter: public OutputFilePack {

	static int ecode(const Element::CODE &code)
	{
		switch (code) {
		case Element::CODE::POINT1:
			return 1;
		case Element::CODE::LINE2:
			return 3;
		case Element::CODE::LINE3:
			return 21;
		case Element::CODE::SQUARE4:
			return 9;
		case Element::CODE::SQUARE8:
			return 23;
		case Element::CODE::TRIANGLE3:
			return 5;
		case Element::CODE::TRIANGLE6:
			return 22;
		case Element::CODE::TETRA4:
			return 10;
		case Element::CODE::TETRA10:
			return 24;
		case Element::CODE::PYRAMID5:
			return 14;
		case Element::CODE::PYRAMID13:
			return 27;
		case Element::CODE::PRISMA6:
			return 13;
		case Element::CODE::PRISMA15:
			return 26;
		case Element::CODE::HEXA8:
			return 12;
		case Element::CODE::HEXA20:
			return 25;
		default:
			return -1;
		}
	}

	void description(const std::string &description)
	{
		insert(snprintf(buffer, bsize, "%s", description.c_str()));
	}

	void push(char c)
	{
		insert(1, c);
	}

	void points(int size)
	{
		insert(snprintf(buffer, bsize, "\nPOINTS %d float\n", size));
	}

	void point(float x, float y, float z)
	{
		insert(snprintf(buffer, bsize, "%12.5e %12.5e %12.5e\n", x, y, z));
	}

	void cells(int cells, int size)
	{
		insert(snprintf(buffer, bsize, "\nCELLS %d %d\n", cells, size));
	}

	void cell(int size, const esint* nodes)
	{
		insert(snprintf(buffer, bsize, "%d", size));
		for (int n = 0; n < size; ++n) {
			insert(snprintf(buffer, bsize, " %d", (int)nodes[n]));
		}
		insert(snprintf(buffer, bsize, "\n"));
	}

	void cell(int size, const esint* nodes, esint offset)
	{
		insert(snprintf(buffer, bsize, "%d", size));
		for (int n = 0; n < size; ++n) {
			insert(snprintf(buffer, bsize, " %d", (int)(nodes[n] + offset)));
		}
		insert(snprintf(buffer, bsize, "\n"));
	}

	void cell(int size, const esint* nodes, const esint* map)
	{
		insert(snprintf(buffer, bsize, "%d", size));
		for (int n = 0; n < size; ++n) {
			insert(snprintf(buffer, bsize, " %d", (int)map[nodes[n]]));
		}
		insert(snprintf(buffer, bsize, "\n"));
	}

	void cell(int size, const esint* nodes, const esint* map,  esint offset)
	{
		insert(snprintf(buffer, bsize, "%d", size));
		for (int n = 0; n < size; ++n) {
			insert(snprintf(buffer, bsize, " %d", (int)map[nodes[n]] + offset));
		}
		insert(snprintf(buffer, bsize, "\n"));
	}

	void type(const Element::CODE &code)
	{
		insert(snprintf(buffer, bsize, "%d\n", ecode(code)));
	}

	void celltypes(int size)
	{
		insert(snprintf(buffer, bsize, "\nCELL_TYPES %d\n", size));
	}

	void pointdata(int size)
	{
		insert(snprintf(buffer, bsize, "\nPOINT_DATA %d\n", size));
	}

	void celldata(int size)
	{
		insert(snprintf(buffer, bsize, "\nCELL_DATA %d\n", size));
	}

	void data(const std::string &type, const std::string &name, const std::string &format)
	{
		insert(snprintf(buffer, bsize, "%s %s %s\n", type.c_str(), name.c_str(), format.c_str()));
	}

	void int32s(int value)
	{
		insert(snprintf(buffer, bsize, "%d ", value));
	}

	void int32ln(int value)
	{
		insert(snprintf(buffer, bsize, "%d\n", value));
	}

	void float32s(float value)
	{
		insert(snprintf(buffer, bsize, "%f ", value));
	}

	void float32ln(float value)
	{
		insert(snprintf(buffer, bsize, "%f\n", value));
	}
};
}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_VTKWRITTER_H_ */
