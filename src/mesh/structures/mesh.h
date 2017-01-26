#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>
#include <functional>

#include "coordinates.h"
#include "material.h"
#include "region.h"

#include "esbasis.h"

namespace espreso {

struct Configuration;
template <typename TParameter, typename TValue>
struct ConfigurationVectorMap;

class Element;

namespace input {
class Loader;
}

class Mesh
{

public:
	friend class input::Loader;

	Mesh();
	virtual ~Mesh();

	virtual void partitiate(size_t parts);
	void computeFixPoints(size_t number);
	void computeVolumeCorners(size_t number, bool onVertices, bool onEdges, bool onFaces);
	void computePlaneCorners(size_t number, bool onVertices, bool onEdges);

	void computeElementsFromFaces();

	void computeFacesOfAllElements();
	void computeFacesOnDomainsSurface();
	virtual void computeFacesSharedByDomains();
	void clearFacesWithoutSettings();

	void computeEdgesOfAllElements();
	void computeEdgesSharedByDomains();
	void computeEdgesOnBordersOfFacesSharedByDomains();
	void clearEdgesWithoutSettings();

	void computeCornersOnEdges(size_t number, bool onVertices, bool onEdges);
	void computeCornersOnFaces(size_t number, bool onVertices, bool onEdges, bool onFaces);

	void loadProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep = 0);
	void loadNodeProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep = 0);
	void loadProperty(const ConfigurationVectorMap<std::string, std::string> &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties);
	void loadNodeProperty(const ConfigurationVectorMap<std::string, std::string> &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties);

	void removeDuplicateRegions();

	template<typename TMaterial>
	void loadMaterials(const std::map<std::string, TMaterial*> &materials, const std::map<std::string, std::string> &sets)
	{
		size_t index = 0;
		for (auto it = sets.begin(); it != sets.end(); ++it, index++) {
			Region *region = this->region(it->first);
			#pragma omp parallel for
			for (size_t e = 0; e < region->elements().size(); e++) {
				region->elements()[e]->setParam(Element::MATERIAL, index);
			}
			_materials.push_back(new Material(_coordinates, *materials.find(it->second)->second));
			ESINFO(OVERVIEW) << "Set material '" << it-> second << "' for region '" << region->name << "'";
		}
		if (!_materials.size()) {
			ESINFO(GLOBAL_ERROR) << "ESPRESO needs at least one material.";
		}
	}

	const Coordinates& coordinates() const { return _coordinates; }
	const std::vector<Element*>& elements() const { return _elements; };
	const std::vector<Element*>& faces() const { return _faces; };
	const std::vector<Element*>& edges() const { return _edges; };
	const std::vector<Element*>& nodes() const { return _nodes; };

	const Element* getDOFsElement(size_t part, size_t DOF) const { return _DOFtoElement[part][DOF]; }

	const std::vector<Element*>& fixPoints(size_t part) const { return _fixPoints[part]; };
	const std::vector<Element*>& corners() const { return _corners; };

	size_t parts() const { return _partPtrs.size() - 1; }
	const std::vector<eslocal>& getPartition() const { return _partPtrs; }

	const std::vector<int>& neighbours() const { return _neighbours; }
	const std::vector<Region*>& regions() const { return _regions; }
	const std::vector<Material*>& materials() const { return _materials; }
	const std::vector<Evaluator*>& evaluators() const { return _evaluators; }

	std::vector<size_t> assignVariousDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToEdges(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToFaces(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToElements(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);

	void computeNodesDOFsCounters(const std::vector<Property> &DOFs);
	void computeEdgesDOFsCounters(const std::vector<Property> &DOFs);
	void computeFacesDOFsCounters(const std::vector<Property> &DOFs);

	void getSurface(Mesh &surface) const;

	void synchronizeGlobalIndices();

	Region* region(const std::string &name)
	{
		auto it = std::find_if(_regions.begin(), _regions.end(), [&] (const Region *region) { return region->name.compare(name) == 0; });
		if (it != _regions.end()) {
			return *it;
		}
		ESINFO(GLOBAL_ERROR) << "Unknown region '" << name << "'";
		exit(EXIT_FAILURE);
	}

protected:
	void fillFacesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* face)> filter);
	void fillEdgesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* edge)> filter);
	void fillNodesFromCoordinates();

	void fillEdgesFromFaces(std::function<bool(const std::vector<Element*> &nodes, const Element* edge)> filter);

	void fillParentEdgesToNodes();
	void fillParentFacesToNodes();
	void fillParentElementsToNodes();

	void fillEdgesParents();
	void fillFacesParents();

	void mapFacesToClusters();
	void mapEdgesToClusters();

	void mapCoordinatesToDomains();
	void mapElementsToDomains();
	void mapFacesToDomains();
	void mapEdgesToDomains();
	void mapNodesToDomains();

	std::vector<eslocal> getPartition(size_t begin, size_t end, eslocal parts) const;
	std::vector<eslocal> getPartition(const std::vector<Element*> &elements, size_t begin, size_t end, eslocal parts) const;
	eslocal getCentralNode(eslocal begin, eslocal end, const std::vector<eslocal> &ePartition, eslocal part, eslocal subpart) const;
	eslocal getCentralNode(const std::vector<Element*> &elements, size_t begin, size_t end, const std::vector<eslocal> &ePartition, eslocal subpart) const;
	void makePartContinuous(size_t part);

	/** @brief Reference to coordinates. */
	Coordinates _coordinates;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<eslocal> _partPtrs;

	/// Elements of the mesh.
	std::vector<Element*> _elements;

	/// Faces of the elements.
	std::vector<Element*> _faces;

	/// Edges of the elements.
	std::vector<Element*> _edges;

	/// Nodes of the elements.
	std::vector<Element*> _nodes;

	/// Back map from DOFs to Element
	std::vector<std::vector<Element*> > _DOFtoElement;

	/** @brief Fix points for all parts. */
	std::vector<std::vector<Element*> > _fixPoints;

	/// Corners for HFETI
	std::vector<Element*> _corners;

	/** @brief list of neighbours MPI ranks */
	std::vector<int> _neighbours;

	/** @brief list of materials in the mesh*/
	std::vector<Material*> _materials;

	/** @brief list of mesh regions*/
	std::vector<Region*> _regions;

	/** @brief list of evaluators */
	std::vector<Evaluator*> _evaluators;

private:
	Mesh(const Mesh &mesh)
	{
		ESINFO(ERROR) << "It is not allowed to copy Mesh.";
	}

	Mesh& operator=(const Mesh &mesh)
	{
		ESINFO(ERROR) << "It is not allowed to copy Mesh.";
		return *this;
	}

};

namespace input { class API; }

class APIMesh: public Mesh
{
	friend class input::API;
public:
	APIMesh(eslocal *l2g, size_t size): _l2g(l2g, l2g + size) { };

	void partitiate(size_t parts);

	void mapDOFsToDomains();
	void fillParentElementsToDOFs(const std::vector<std::vector<eslocal> > &eDOFs);

	void computeFacesSharedByDomains();
	void computeDOFsDOFsCounters();
	std::vector<size_t> distributeDOFsToDomains(const std::vector<size_t> &offsets);

	const std::vector<Element*>& DOFs() const { return _DOFs; }

protected:
	std::vector<Element*> _DOFs;
	std::vector<esglobal> _l2g;
	std::vector<G2L> _g2l;
};

}


#endif /* MESH_H_ */
