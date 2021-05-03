
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include "element.h"

#include <string>
#include <vector>

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

	static bool convertDatabase();

	static void init();
	static void load();
	static void finish();
	static Element edata[(int)Element::CODE::SIZE];

	Mesh();
	~Mesh();
	void preprocess();
	void preprocessForGUI()
	{
		_withGUI = true;
		preprocess();
	}
	void duplicate();
	void toBuffer();
	void setMaterials();
	void printMeshStatistics();
	void printDecompositionStatistics();

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

	ElementStore *halo;

	SurfaceStore *surface;
	SurfaceStore *domainsSurface;

	ContactStore *contact;

	std::vector<int> neighbors;
	std::vector<int> neighborsWithMe;

	std::vector<const MaterialConfiguration*> materials;

	Output *output;

protected:
	void analyze();
	void computePersistentParameters();

	bool _omitClusterization, _omitDecomposition;
	bool _withGUI, _withFETI, _withBEM, _withEdgeDual;
};

}



#endif /* SRC_MESH_MESH_H_ */
