
#ifndef SRC_INPUT_MESHDATA_H_
#define SRC_INPUT_MESHDATA_H_

#include "basis/containers/point.h"

#include <vector>
#include <utility>
#include <map>

namespace espreso {

// Mesh data for building the mesh in ESPRESO

struct InputTransformationConfiguration;

struct MeshData {
	enum class TYPE {
		GENERAL,
		SORTED, // not implemented yet
		GENERATED,
	};

	MeshData(TYPE type = TYPE::GENERAL): type(type) {}

	TYPE type;
	bool removeDuplicates = false;

	std::vector<esint> nIDs;       // nodes IDs [arbitrary numbers]
	std::vector<Point> coordinates;  // nodes coordinates

	std::vector<esint> eIDs;       // elements IDs [arbitrary numbers]
	std::vector<esint> esize;      // the number of nodes for a given element [4, 4, 4]         e0 has 4 nodes
	std::vector<esint> enodes;     // elements nodes                          [0, 1, 5, 4, ...] e0 has nodes 0, 1, 5, 4
	std::vector<int> etype;          // elements types [from Element::CODE]
	std::vector<int> body;           // elements bodies
	std::vector<int> material;       // elements materials

	std::map<std::string, std::vector<esint> > eregions; // elements regions <name, list of IDs>
	std::map<std::string, std::vector<esint> > nregions; // nodes regions <name, list of IDs>
};

struct MeshBuilder: public MeshData {

	struct Duplicate {
		esint id, target;
		esint idoffset, toffset;

		Duplicate(): id(0), target(0), idoffset(0), toffset(0) {}
		Duplicate(esint id, esint target): id(id), target(target), idoffset(0), toffset(0) {}
		Duplicate(esint id, esint target, esint idoffset, esint toffset): id(id), target(target), idoffset(idoffset), toffset(toffset) {}

		bool operator()(Duplicate const &a, Duplicate const &b) { return a.id < b.id; };
	};

	struct Geometry { // for transformations
		esint nids, eids;
		std::pair<size_t, size_t> nodes, elements, enodes;
		std::vector<std::pair<size_t, size_t>> nregsize, eregsize;
	};

	MeshBuilder(TYPE type = TYPE::GENERAL): MeshData(type) {}
	virtual ~MeshBuilder() {}

	virtual void load() = 0;
	virtual void build();

	void selectGeometry(Geometry &geometry);
	void translate(InputTransformationConfiguration &transformation, Geometry &source, int &instance);
	void rotate(InputTransformationConfiguration &transformation, Geometry &source, int &instance);
	void scale(InputTransformationConfiguration &transformation, Geometry &source, int &instance);
	void shear(InputTransformationConfiguration &transformation, Geometry &source, int &instance);

	std::vector<esint> _nrankdist; // nodes ranks distribution [0, 2, 5, ...] n0 is on 2 processes
	std::vector<int> _nranks;        // nodes ranks              [0, 1, ...]    n0 is on processes 0 and 1
	std::vector<esint> _edist;     // elements nodes distribution [0, 4, 8, ...]

	std::vector<Duplicate> _duplicateNodes;
	std::vector<Duplicate> _duplicateElements;

private:
	void duplicate(Geometry &source, int instance);
};

}



#endif /* SRC_INPUT_MESHDATA_H_ */

