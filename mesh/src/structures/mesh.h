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
#include "esbem.h"

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

	Mesh();

	Mesh(const Mesh &other);
	Mesh& operator=(const Mesh &other);

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

	size_t getFixPointsCount() const
	{
		return _fixPoints.size() / (_partPtrs.size() - 1);
	}

	const std::vector<eslocal>& getFixPoints() const
	{
		return _fixPoints;
	}

	eslocal getPartNodesCount(eslocal part) const
	{
		return _coordinates.localSize(part);
	}

	void elasticity(SparseCSRMatrix<eslocal> &K, SparseCSRMatrix<eslocal> &M, std::vector<double> &f, eslocal part) const
	{
		SparseVVPMatrix<eslocal> _K;
		SparseVVPMatrix<eslocal> _M;
		_elasticity(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}
	void elasticity(SparseCSRMatrix<eslocal> &K, std::vector<double> &f, eslocal part) const
	{
		SparseVVPMatrix<eslocal> _K;
		SparseVVPMatrix<eslocal> _M;
		_elasticity(_K, _M, f, part, false);
		K = _K;
	}
	void heat(SparseCSRMatrix<eslocal> &K, SparseCSRMatrix<eslocal> &M, std::vector<double> &f, eslocal part) const
	{
		SparseVVPMatrix<eslocal> _K;
		SparseVVPMatrix<eslocal> _M;
		_heat(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}

protected:
	static void assign(Mesh &m1, Mesh &m2);

	void _elasticity(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f, eslocal part, bool dynamic) const;
	void _assembleElesticity(
		const Element *e,
		size_t part,
		DenseMatrix &Ke,
		DenseMatrix &Me,
		std::vector<double> &fe,
		std::vector<double> &inertia,
		DenseMatrix &C,
		bool dynamic) const;

	void _integrateElasticity(
		const Element *e,
		SparseVVPMatrix<eslocal> &K,
		SparseVVPMatrix<eslocal> &M,
		std::vector<double> &f,
		const DenseMatrix &Ke,
		const DenseMatrix &Me,
		const std::vector<double> &fe,
		bool dynamic
	) const;

	void _heat(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f, eslocal part, bool dynamic) const;
	void _assembleHeat(
		const Element *e,
		size_t part,
		DenseMatrix &Ke,
		DenseMatrix &Me,
		std::vector<double> &fe,
		DenseMatrix &C,
		bool dynamic) const;

	void _integrateHeat(
		const Element *e,
		SparseVVPMatrix<eslocal> &K,
		SparseVVPMatrix<eslocal> &M,
		std::vector<double> &f,
		const DenseMatrix &Ke,
		const DenseMatrix &Me,
		const std::vector<double> &fe,
		bool dynamic
	) const;

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
	std::vector<eslocal> _fixPoints;

	/** @brief Map of points to sub-domains. */
	Boundaries _subdomainBoundaries;

	/** @brief Map of points to clusters. */
	Boundaries _clusterBoundaries;
};


class SurfaceMesh: public Mesh
{

public:
	SurfaceMesh(): Mesh() { };
	SurfaceMesh(const Mesh &mesh)
	{
		mesh.getSurface(*this);
	}

	void elasticity(DenseMatrix &K, size_t part) const;
	void integrateUpperFaces(std::vector < double > &f , size_t part);
};

class CommonFacesMesh: public Mesh
{

public:
	CommonFacesMesh(): Mesh() { };
	CommonFacesMesh(const Mesh &mesh)
	{
		mesh.getCommonFaces(*this);
	}
};

class CornerLinesMesh: public Mesh
{

public:
	CornerLinesMesh(): Mesh() { };
	CornerLinesMesh(const Mesh &mesh)
	{
		mesh.getCornerLines(*this);
	}
};

}


#endif /* MESH_H_ */
