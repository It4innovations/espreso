
#ifndef SRC_OUTPUT_REGIONDATA_H_
#define SRC_OUTPUT_REGIONDATA_H_

#include <vector>
#include <map>
#include <string>

namespace espreso {

struct Solution;

namespace output {

struct DataArrays {
	std::map<std::string, std::pair<size_t, std::vector<eslocal>* > > pointDataInteger, elementDataInteger;
	std::map<std::string, std::pair<size_t, std::vector<double>* > > pointDataDouble, elementDataDouble;

	void clear();
	~DataArrays();
};

struct RegionData {
	std::vector<double> coordinates; // x1, y1, z1, x2, y2, z2, ...
	std::vector<eslocal> elementsTypes;  // code1, code2, ...
	std::vector<eslocal> elementsNodes;  // nodes1, nodes2, ...
	std::vector<eslocal> elements;  // n11, n12, n13, ..., n21, n22, n23, ...

	DataArrays data;

	std::vector<std::string> pointDataNames() const;
	std::vector<std::string> cellDataNames() const;

	void clearData();

	size_t packedSize() const;
	void pack(char* data) const;
	void pack(std::vector<char> &data) const;
	void unpack(const char* &data);
};

}
}



#endif /* SRC_OUTPUT_REGIONDATA_H_ */
