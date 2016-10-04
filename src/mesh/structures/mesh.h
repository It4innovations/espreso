#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>
#include <functional>

#include "metis.h"

#include "../elements/elements.h"
#include "coordinates.h"
#include "material.h"

#include "esbasis.h"

namespace espreso {

namespace input {
class Loader;
}

class Boundaries;

class Mesh
{

public:
	friend class input::Loader;

	Mesh();
	virtual ~Mesh();

	virtual void partitiate(size_t parts);
	void computeFixPoints(size_t number);
	void computeVolumeCorners(size_t number, bool vertices, bool edges, bool faces);
	void computePlaneCorners(size_t number, bool vertices, bool edges);

	void computeFacesOfAllElements();
	void computeFacesOnDomainsSurface();
	void computeFacesSharedByDomains();
	void clearFacesWithoutSettings();

	void computeEdgesOfAllElements();
	void computeEdgesSharedByDomains();
	void computeEdgesOnBordersOfFacesSharedByDomains();
	void clearEdgesWithoutSettings();

	void computeCornersOnEdges(size_t number);
	void computeCornersOnFaces(size_t number);

	const Coordinates& coordinates() const { return _coordinates; }
	const std::vector<Element*>& elements() const { return _elements; };
	const std::vector<Element*>& faces() const { return _faces; };
	const std::vector<Element*>& edges() const { return _edges; };
	const std::vector<Element*>& nodes() const { return _nodes; };

	const std::vector<Element*>& fixPoints(size_t part) const { return _fixPoints[part]; };
	const std::vector<Element*>& corners() const { return _corners; };

	size_t parts() const { return _partPtrs.size() - 1; }
	const std::vector<eslocal>& getPartition() const { return _partPtrs; }

	const std::vector<int>& neighbours() const { return _neighbours; }
	const std::vector<Material>& materials() const { return _materials; }
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

	/** @brief Fix points for all parts. */
	std::vector<std::vector<Element*> > _fixPoints;

	/// Corners for HFETI
	std::vector<Element*> _corners;

	/** @brief list of neighbours MPI ranks */
	std::vector<int> _neighbours;

	/** @brief list of materials in the mesh*/
	std::vector<Material> _materials;

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
	APIMesh(const std::vector<std::vector<eslocal> > &eDOFs, const std::vector<std::vector<double> > &eMatrices)
	: _eDOFs(eDOFs), _eMatrices(eMatrices) { };

	void partitiate(size_t parts);

	void mapDOFsToDomains();
	void fillParentElementsToDOFs(const std::vector<std::vector<eslocal> > &eDOFs);

	void computeDOFsDOFsCounters();
	std::vector<size_t> distributeDOFsToDomains(const std::vector<size_t> &offsets);

	const std::vector<Element*> DOFs() const { return _DOFs; }

	const std::vector<eslocal>& eDOFs(size_t index) const { return _eDOFs[index]; }
	const std::vector<double>& eMatrix(size_t index) const { return _eMatrices[index]; }
	size_t elementSize() const { return _eMatrices.size(); }

protected:
	std::vector<Element*> _DOFs;
	const std::vector<std::vector<eslocal>> &_eDOFs;
	const std::vector<std::vector<double> > &_eMatrices;
};

}


#endif /* MESH_H_ */
