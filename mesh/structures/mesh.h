#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "mkl_spblas.h"
#include "mkl_blas.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

#include "cilk/cilk.h"

#include "metis.h"

#include "../elements/elements.h"
#include "coordinates.h"
#include "boundaries.h"

#include "esbasis.h"

namespace esinput {
class InternalLoader;
class ExternalLoader;
}

namespace mesh {


class Boundaries;

class SurfaceMesh;
class CommonFacesMesh;
class CornerLinesMesh;

class Mesh
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);
	friend class esinput::InternalLoader;
	friend class esinput::ExternalLoader;

	Mesh(int rank, int size);

	const Coordinates& coordinates() const
	{
		return _coordinates;
	}

	Coordinates& coordinates()
	{
		return _coordinates;
	}

	const Boundaries& subdomainBoundaries() const
	{
		return _subdomainBoundaries;
	}

	const Boundaries& clusterBoundaries() const
	{
		return _clusterBoundaries;
	}

	int rank() const
	{
		return _rank;
	}

	int size() const
	{
		return _size;
	}

	void saveNodeArray(eslocal *nodeArray, size_t part);

	void getSurface(SurfaceMesh &surface) const;
	void getCommonFaces(CommonFacesMesh &commonFaces) const;
	void getCornerLines(CornerLinesMesh &cornerLines) const;

	~Mesh();

	void partitiate(eslocal parts, eslocal fixPoints);
	void computeFixPoints(eslocal fixPoints);
	void computeCorners(Boundaries &boundaries, eslocal number, bool corners, bool edges, bool faces) const;

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

protected:
	eslocal* getPartition(eslocal first, eslocal last, eslocal parts) const;
	eslocal getCentralNode(eslocal first, eslocal last, eslocal *ePartition, eslocal part, eslocal subpart) const;

	void partitiate(eslocal *ePartition);
	void computeLocalIndices(size_t part);
	void computeBoundaries();

	void checkMETISResult(eslocal result) const;
	void checkMKLResult(eslocal result) const;

	/** @brief Reference to coordinates. */
	Coordinates _coordinates;

	/** @brief Array that stores all elements of the mesh. */
	std::vector<mesh::Element*> _elements;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<eslocal> _partPtrs;

	/** @brief Fix points for all parts. */
	std::vector<std::vector<eslocal> > _fixPoints;

	/** @brief Map of points to sub-domains. */
	Boundaries _subdomainBoundaries;

	/** @brief Map of points to clusters. */
	Boundaries _clusterBoundaries;

	/** @brief MPI rank */
	int _rank;

	/** @brief MPI size */
	int _size;
};


class SurfaceMesh: public Mesh
{

public:
	SurfaceMesh(int rank, int size): Mesh(rank, size) { };
	SurfaceMesh(const Mesh &mesh): Mesh(mesh.rank(), mesh.size())
	{
		mesh.getSurface(*this);
	}
};

class CommonFacesMesh: public Mesh
{

public:
	CommonFacesMesh(int rank, int size): Mesh(rank, size) { };
	CommonFacesMesh(const Mesh &mesh): Mesh(mesh.rank(), mesh.size())
	{
		mesh.getCommonFaces(*this);
	}
};

class CornerLinesMesh: public Mesh
{

public:
	CornerLinesMesh(int rank, int size): Mesh(rank, size) { };
	CornerLinesMesh(const Mesh &mesh): Mesh(mesh.rank(), mesh.size())
	{
		mesh.getCornerLines(*this);
	}
};

}


#endif /* MESH_H_ */
