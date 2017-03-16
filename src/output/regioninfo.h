
#ifndef SRC_OUTPUT_REGIONINFO_H_
#define SRC_OUTPUT_REGIONINFO_H_

#include <vector>
#include <map>
#include <string>

namespace espreso {

class Mesh;
class Element;
struct Region;
struct Solution;
struct Point;

namespace output {

struct DataArrays {
	std::map<std::string, std::pair<size_t, std::vector<eslocal>* > > pointDataInteger, elementDataInteger;
	std::map<std::string, std::pair<size_t, std::vector<double>* > > pointDataDouble, elementDataDouble;

	void clear();
	~DataArrays();
};

struct RegionInfo {

	RegionInfo(const Mesh *mesh): _mesh(mesh), _body(-1), _region(NULL) {};
	RegionInfo(const Mesh *mesh, size_t body): _mesh(mesh), _body(body), _region(NULL) {};
	RegionInfo(const Mesh *mesh, const Region *region): _mesh(mesh), _body(-1), _region(region) {};

	virtual RegionInfo* deriveRegion(const Region *region) const =0;
	virtual RegionInfo* copyWithoutMesh() const =0;

	virtual ~RegionInfo() {};

	virtual void addSettings(size_t step) =0;
	virtual void addSolution(const std::vector<Solution*> &solution) =0;
	virtual void addGeneralInfo() = 0;

	void clearData();

	virtual bool isShrunk() const =0;
	virtual bool distributed() const =0;
	virtual Point shrink(const Point &p, eslocal domain) const =0;

	std::vector<double> coordinates; // x1, y1, z1, x2, y2, z2, ...
	std::vector<eslocal> elementsTypes;  // code1, code2, ...
	std::vector<eslocal> elementsNodes;  // nodes1, nodes2, ...
	std::vector<eslocal> elements;  // n11, n12, n13, ..., n21, n22, n23, ...

	DataArrays data;
	std::vector<Solution*> solutions;

protected:
	const Mesh *_mesh;
	size_t _body;
	const Region *_region;
};

}
}




#endif /* SRC_OUTPUT_REGIONINFO_H_ */
