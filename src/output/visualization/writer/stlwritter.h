
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_

#include "basis/io/outputfile.h"

namespace espreso {

struct STLBinaryOutputWriter: public OutputFilePack {

	void storeHeader(const std::string &value)
	{
		insert(value.size(), value.data());
		insert(80 - value.size(), '\0');
	}

	void storeFooter(const std::string &value)
	{

	}

	void storeSize(int size)
	{
		insert(sizeof(int), &size);
	}

	void beginFace(float x, float y, float z)
	{
		insert(sizeof(float), &x);
		insert(sizeof(float), &y);
		insert(sizeof(float), &z);
	}

	void addVertex(float x, float y, float z)
	{
		insert(sizeof(float), &x);
		insert(sizeof(float), &y);
		insert(sizeof(float), &z);
	}

	void endFace()
	{
		int attribute = 0;
		insert(2, &attribute);
	}
};

struct STLASCIIOutputWriter: public OutputFilePack {

	inline void storeHeader(const std::string &name)
	{
		insert(snprintf(buffer, bsize, "solid %s\n", name.c_str()));
	}

	inline void storeFooter(const std::string &name)
	{
		insert(snprintf(buffer, bsize, "endsolid %s\n", name.c_str()));
	}

	inline void storeSize(int size)
	{

	}

	inline void beginFace(float x, float y, float z)
	{
		insert(snprintf(buffer, bsize, " facet normal %f %f %f\n  outer loop\n", x, y ,z));
	}

	inline void addVertex(float x, float y, float z)
	{
		insert(snprintf(buffer, bsize, "   vertex %f %f %f\n", x, y ,z));
	}

	inline void endFace()
	{
		insert(snprintf(buffer, bsize, "  endloop\n endfacet\n"));
	}
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_ */
