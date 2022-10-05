
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include "element.h"
#include "basis/containers/point.h"

#include <string>
#include <vector>
#include <map>
#include <mutex>

namespace espreso {

struct ECF;
struct OutputConfiguration;
struct MaterialConfiguration;

struct Statistics;
struct ElementStore;
struct ElementData;
struct NodeStore;
struct DomainStore;
struct ClusterStore;
struct BodyStore;
struct NodeData;

struct ElementsRegionStore;
struct BoundaryRegionStore;
struct ContactInterfaceStore;
struct FETIDataStore;
struct SurfaceStore;
struct ContactStore;

class Output;

class Mesh {
public:
	static void init();
	static void load();
	static void finish();
	static Element edata[(int)Element::CODE::SIZE];
	template <typename TCode>
	static inline const Element& element(const TCode &code)
	{
		return edata[(short)code & 255];
	}

	Mesh();
	~Mesh();
	void preprocess();
	void preprocessForGUI()
	{
		_withGUI = true;
		preprocess();
	}
	void updateMesh();
	void partitiate(int ndomains);
	void duplicate();
	void toBuffer();
	void printMeshStatistics();
	void printDecompositionStatistics();
	size_t meshSize();

	ElementsRegionStore* allElements()
	{
		return eregion("ALL_ELEMENTS");
	}

	BoundaryRegionStore* allNodes()
	{
		return bregion("ALL_NODES");
	}

	ElementsRegionStore* eregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);

	bool onAllElements(const std::string &eregion) const;

	bool hasPhaseChange() const;

	int dimension;
	size_t preferedDomains;

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;
	std::vector<ContactInterfaceStore*> contactInterfaces;

	DomainStore *domains;
	ClusterStore *clusters;
	BodyStore *bodies;

	FETIDataStore *FETIData;

	SurfaceStore *surface;
	SurfaceStore *domainsSurface;

	ContactStore *contact;

	std::vector<int> neighbors;
	std::vector<int> neighborsWithMe;

	std::vector<const MaterialConfiguration*> materials;
	std::map<std::string, _Point<esfloat> > orientation;

	Output *output;

	bool withDualGraph, withSeparatedRegions, convertToVolume, balanceVolume;

protected:
	void analyze();
	void setMaterials();
	void reclusterize();
	void computePersistentParameters();

	bool _omitClusterization, _omitDecomposition, _omitDual, _omitRegionInfo, _omitSynchronization;
	bool _withGUI, _withFETI, _withBEM, _withEdgeDual;
};

}



#endif /* SRC_MESH_MESH_H_ */
