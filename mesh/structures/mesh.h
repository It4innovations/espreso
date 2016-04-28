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

#include "metis.h"

#include "../elements/elements.h"
#include "coordinates.h"
#include "boundaries.h"
#include "faces.h"
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

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);
	friend class input::Loader;

	Mesh();

	const Coordinates& coordinates() const
	{
		return _coordinates;
	}

	Coordinates& coordinates()
	{
		return _coordinates;
	}

	Faces& faces() {
		return _faces;
	}

	const Boundaries& subdomainBoundaries() const
	{
		return _subdomainBoundaries;
	}

	const Boundaries& clusterBoundaries() const
	{
		return _clusterBoundaries;
	}

	void saveNodeArray(eslocal *nodeArray, size_t part) const;

	void getSurface(Mesh &surface) const;

	virtual ~Mesh();

	virtual void partitiate(size_t parts);
	void computeFixPoints(size_t fixPoints);
	void computeCorners(eslocal number, bool vertices, bool edges, bool faces, bool averageEdges, bool averageFaces);

	const std::vector<Element*>& getElements() const
	{
		return _elements;
	};

	size_t parts() const
	{
		return _partPtrs.size() - 1;
	}

	const std::vector<eslocal>& getPartition() const
	{
		return _partPtrs;
	}

	const std::vector<std::vector<eslocal> >& getFixPoints() const
	{
		return _fixPoints;
	}

	eslocal getPartNodesCount(eslocal part) const
	{
		return _coordinates.localSize(part);
	}

	const std::vector<int>& neighbours() const
	{
		return _neighbours;
	}

	size_t DOFs() const
	{
		return _DOFs;
	}

	const std::vector<Material>& materials() const
	{
		return _materials;
	}

protected:
	eslocal* getPartition(eslocal first, eslocal last, eslocal parts) const;
	eslocal getCentralNode(eslocal first, eslocal last, eslocal *ePartition, eslocal part, eslocal subpart) const;

	void computeBoundaries();
	void remapElementsToSubdomain();
	void remapElementsToCluster();

	void makePartContinuous(size_t part);
	void computeCommonFaces(Mesh &faces);
	void computeBorderLinesAndVertices(const Mesh &faces, std::vector<bool> &border, Mesh &lines, std::set<eslocal> &vertices);
	void prepareAveragingLines(Mesh &faces, Mesh &lines);
	void prepareAveragingFaces(Mesh &faces, std::vector<bool> &border);
	void correctCycle(Mesh &faces, Mesh &lines, bool average);

	/** @brief Reference to coordinates. */
	Coordinates _coordinates;

	/** @brief Indexes of all faces, it's element contains an element reference.*/
	Faces _faces;

	/** @brief Array that stores all elements of the mesh. */
	std::vector<Element*> _elements;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<eslocal> _partPtrs;

	/** @brief Fix points for all parts. */
	std::vector<std::vector<eslocal> > _fixPoints;

	/** @brief Map of points to sub-domains. */
	Boundaries _subdomainBoundaries;

	/** @brief Map of points to clusters. */
	Boundaries _clusterBoundaries;

	/** @brief list of neighbours MPI ranks */
	std::vector<int> _neighbours;

	/** @brief the number of DOFs for all nodes*/
	size_t _DOFs;

	/** @brief list of materials in the mesh*/
	std::vector<Material> _materials;

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
