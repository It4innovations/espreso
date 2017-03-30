
#ifndef SRC_OUTPUT_MESHINFO_H_
#define SRC_OUTPUT_MESHINFO_H_

#include <vector>
#include <map>
#include <string>

#include "regiondata.h"

namespace espreso {

class Mesh;
struct Region;
struct Solution;
struct Point;

namespace output {

class MeshInfo {

	friend class ResultStore;

public:

	enum InfoMode: int {
		EMPTY              = 0,
		PREPARE            = 1 << 0,
		SEPARATE_BODIES    = 1 << 1,
		SEPARATE_MATERIALS = 1 << 2
	};

	MeshInfo(const Mesh *mesh, InfoMode mode): _mesh(mesh), _mode(mode), _region(NULL) {};
	MeshInfo(const Mesh *mesh, const Region *region, InfoMode mode): _mesh(mesh), _mode(mode), _region(region) {};

	size_t regions() const { return _regions.size(); }
	const RegionData& region(size_t r) const { return _regions[r]; }

	virtual MeshInfo* deriveRegion(const Region *region) const =0;
	virtual MeshInfo* copyWithoutMesh() const =0;

	virtual ~MeshInfo() {};

	virtual void addSettings(size_t step) =0;
	virtual void addSolution(const std::vector<Solution*> &solution) =0;
	virtual void addGeneralInfo() = 0;

	void clearData();

	virtual bool isShrunk() const =0;
	virtual bool distributed() const =0;
	virtual Point shrink(const Point &p, eslocal domain) const =0;

protected:
	const Mesh *_mesh;
	InfoMode _mode;

	std::vector<RegionData> _regions;
	const Region *_region;
};

inline MeshInfo::InfoMode operator~(MeshInfo::InfoMode m)
{
	return static_cast<MeshInfo::InfoMode>(~static_cast<int>(m));
}

inline MeshInfo::InfoMode operator|(MeshInfo::InfoMode m1, const MeshInfo::InfoMode &m2)
{
	return static_cast<MeshInfo::InfoMode>(static_cast<int>(m1) | static_cast<int>(m2));
}

inline MeshInfo::InfoMode& operator|=(MeshInfo::InfoMode &m1, const MeshInfo::InfoMode &m2)
{
	m1 = static_cast<MeshInfo::InfoMode>(static_cast<int>(m1) | static_cast<int>(m2));
	return m1;
}

inline MeshInfo::InfoMode operator&(MeshInfo::InfoMode m1, const MeshInfo::InfoMode &m2)
{
	return static_cast<MeshInfo::InfoMode>(static_cast<int>(m1) & static_cast<int>(m2));
}

inline MeshInfo::InfoMode& operator&=(MeshInfo::InfoMode &m1, const MeshInfo::InfoMode &m2)
{
	m1 = static_cast<MeshInfo::InfoMode>(static_cast<int>(m1) & static_cast<int>(m2));
	return m1;
}

}
}




#endif /* SRC_OUTPUT_MESHINFO_H_ */
