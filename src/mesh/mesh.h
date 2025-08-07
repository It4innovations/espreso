
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include "element.h"
#include "basis/containers/point.h"

#include <string>
#include <vector>
#include <map>

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
struct DomainSurfaceStore;
struct ContactStore;

class Output;

class Mesh {
public:

    static bool convertDatabase();

    static void init();
    static void finish();
    static Element edata[(int)Element::CODE::SIZE];

    Mesh();
    ~Mesh();
    void load();
    void preprocess();
    void preprocessForGUI()
    {
        _withGUI = true;
        preprocess();
    }
    void updateMeshCoordinates(const double *displacement);
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

    size_t eregionIndex(const std::string &name);
    size_t bregionIndex(const std::string &name);

    bool onAllElements(const std::string &eregion) const;

    bool hasPhaseChange() const;

    int dimension;
    size_t preferedDomains;

    ElementStore* elements;
    NodeStore* nodes;

    std::vector<ElementsRegionStore*> elementsRegions;
    std::vector<BoundaryRegionStore*> boundaryRegions;
    std::vector<ContactInterfaceStore*> contactInterfaces;
    std::vector<BoundaryRegionStore*> boundary; // boundaryRegions + contactInterfaces

    DomainStore *domains;
    ClusterStore *clusters;
    BodyStore *bodies;

    FETIDataStore *FETIData;

    SurfaceStore *surface;
    DomainSurfaceStore *domainsSurface;

    ContactStore *contact;

    std::vector<int> neighbors;
    std::vector<int> neighborsWithMe;

    std::vector<MaterialConfiguration*> materials;
    std::map<std::string, Point > orientation;

    Output *output;

    bool withSurface;

protected:
    void analyze();
    void setMaterials();
    void reclusterize();
    void computePersistentParameters();

    bool _omitClusterization, _omitDecomposition;
    bool _withGUI, _withFETI, _withBEM, _withEdgeDual;
};

}



#endif /* SRC_MESH_MESH_H_ */
