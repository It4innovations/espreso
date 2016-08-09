#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <functional>

#include "mkl_spblas.h"
#include "mkl_blas.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

#include "cilk/cilk.h"
#include "mkl.h"

#include "metis.h"

#include "../elements/elements.h"
#include "coordinates.h"
#include "material.h"

#include "esbasis.h"
#include "esconfig.h"

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
	void prepare(bool faces, bool edges);
	std::vector<eslocal> computeFixPoints(size_t part, size_t number) const;
	void computeCorners(eslocal number, bool vertices, bool edges, bool faces, bool averageEdges, bool averageFaces);

	const Coordinates& coordinates() const { return _coordinates; }
	const std::vector<Element*>& elements() const { return _elements; };
	const std::vector<Element*>& faces() const { return _faces; };
	const std::vector<Element*>& edges() const { return _edges; };
	const std::vector<Element*>& nodes() const { return _nodes; };

	size_t parts() const { return _partPtrs.size() - 1; }
	const std::vector<eslocal>& getPartition() const { return _partPtrs; }

	const std::vector<int>& neighbours() const { return _neighbours; }
	const std::vector<Material>& materials() const { return _materials; }

	std::vector<size_t> assignVariousDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToEdges(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToFaces(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);
	std::vector<size_t> assignUniformDOFsIndicesToElements(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs);

	void connectNodesDOFsAmongClusters(const std::vector<Property> &DOFs);
	void connectEdgesDOFsAmongClusters(const std::vector<Property> &DOFs);
	void connectFacesDOFsAmongClusters(const std::vector<Property> &DOFs);

	void saveNodeArray(eslocal *nodeArray, size_t part) const;
	void getSurface(Mesh &surface) const;
	std::vector<std::vector<eslocal> > subdomainsInterfaces(Mesh &interface) const;

protected:
	void fillFacesFromElements();
	void fillEdgesFromElements();
	void fillNodesFromElements();

	void mapFacesToClusters();
	void mapEdgesToClusters();

	void mapElementsToDomains();
	void mapFacesToDomains();
	void mapEdgesToDomains();
	void mapNodesToDomains();


	eslocal* getPartition(eslocal first, eslocal last, eslocal parts) const;
	eslocal getCentralNode(eslocal first, eslocal last, eslocal *ePartition, eslocal part, eslocal subpart) const;

	void remapElementsToSubdomain() const;
	void remapElementsToCluster() const;

	void makePartContinuous(size_t part);
	void computeCommonFaces(Mesh &faces);
	void computeBorderLinesAndVertices(const Mesh &faces, std::vector<bool> &border, Mesh &lines, std::set<eslocal> &vertices);
	void prepareAveragingLines(Mesh &faces, Mesh &lines);
	void prepareAveragingFaces(Mesh &faces, std::vector<bool> &border);
	void correctCycle(Mesh &faces, Mesh &lines, bool average);

	/** @brief Reference to coordinates. */
	mutable Coordinates _coordinates;

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
	mutable std::vector<std::vector<eslocal> > _fixPoints;

	/// Corners for HFETI
	std::vector<eslocal> _corners;

	/** @brief list of neighbours MPI ranks */
	std::vector<int> _neighbours;

	/** @brief the number of DOFs for all nodes*/
	size_t _DOFs;

	/** @brief list of materials in the mesh*/
	std::vector<Material> _materials;

	/** @brief list of evaluators */
	std::vector<Evaluator*> _evaluators;

private:
	Mesh(const Mesh &mesh): _DOFs(mesh._DOFs)
	{
		ESINFO(ERROR) << "It is not allowed to copy Mesh.";
	}

	Mesh& operator=(const Mesh &mesh)
	{
		ESINFO(ERROR) << "It is not allowed to copy Mesh.";
		return *this;
	}

};

class APIMesh: public Mesh
{

public:
	APIMesh(std::vector<std::vector<double> > &eMatrices): _eMatrices(eMatrices) { };

	void partitiate(size_t parts);

	std::vector<std::vector<double> >& getMatrices() const
	{
		return _eMatrices;
	}

protected:
	std::vector<std::vector<double> > &_eMatrices;
};

}


#endif /* MESH_H_ */
