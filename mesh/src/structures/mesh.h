#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>

#include "mkl.h"
#include "cilk/cilk.h"

#include "../elements/elements.h"
#include "coordinates.h"

#include "../matrices/sparseDOKMatrix.h"
#include "../settings.h"

#include "esbem.h"

namespace flags
{

enum FLAGS {
	PARTITIONS,
	FIX_POINTS,
	BOUNDARIES,
	NEW_PARTITION,
	FLAGS_SIZE
};
}

class SurfaceMesh;
class CommonFacesMesh;
class CornerLinesMesh;

class Mesh
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);

	Mesh();
	Mesh(const char *mesh, const char *coordinates, idx_t parts, idx_t fixPoints);

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

	void saveVTK(const char* filename, double shrinking = 1);
	void saveVTK(std::vector<std::vector<double> > &displacement, std::vector<std::vector <int> > &l2g_vec );

	void saveNodeArray(idx_t *nodeArray, size_t part);

	void getSurface(SurfaceMesh &surface) const;
	void getCommonFaces(CommonFacesMesh &commonFaces) const;
	void getCornerLines(CornerLinesMesh &cornerLines) const;

	void reserve(size_t size);
	void pushElement(Element* e);
	void endPartition();

	~Mesh();

	void partitiate(idx_t parts, idx_t fixPoints);
	void computeFixPoints(idx_t fixPoints);

	const std::vector<Element*>& getElements() const
	{
		return _elements;
	};

	size_t getPartsCount() const
	{
		return _partsNodesCount.size();
	}

	const std::vector<idx_t>& getPartition() const
	{
		return _partPtrs;
	}

	size_t getFixPointsCount() const
	{
		return _fixPoints.size() / (_partPtrs.size() - 1);
	}

	const std::vector<idx_t>& getFixPoints() const
	{
		return _fixPoints;
	}

	idx_t getPartNodesCount(idx_t part) const
	{
		return _partsNodesCount[part];
	}

	void elasticity(SparseCSRMatrix &K, SparseCSRMatrix &M, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_elasticity(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}
	void elasticity(SparseCSRMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_elasticity(_K, _M, f, part, false);
		K = _K;
	}

protected:
	static void assign(Mesh &m1, Mesh &m2);

	void saveBasis(std::ofstream &vtk, std::vector<std::vector<int> > &l2g_vec);

	Element* createElement(idx_t *indices, idx_t n);

	void _elasticity(SparseVVPMatrix &K, SparseVVPMatrix &M, std::vector<double> &f, idx_t part, bool dynamic);
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
		SparseVVPMatrix &K,
		SparseVVPMatrix &M,
		std::vector<double> &f,
		const DenseMatrix &Ke,
		const DenseMatrix &Me,
		const std::vector<double> &fe,
		bool dynamic
	) const;

	idx_t* getPartition(idx_t first, idx_t last, idx_t parts) const;
	idx_t getCentralNode(idx_t first, idx_t last, idx_t *ePartition, idx_t part, idx_t subpart) const;

	void partitiate(idx_t *ePartition);
	void computeLocalIndices(size_t part);

	void checkMETISResult(int result) const;
	void checkMKLResult(MKL_INT result) const;

	/** @brief Reference to coordinates. */
	Coordinates _coordinates;

	/** @brief The type of indices in element. [local, global] */
	Element::IndicesType _indicesType;

	/** @brief Array that stores all elements of the mesh. */
	std::vector<Element*> _elements;

	/** @brief The biggest node's id in the mesh. */
	idx_t _lastNode;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<idx_t> _partPtrs;

	/** @brief Number of nodes in each part. */
	std::vector<idx_t> _partsNodesCount;

	/** @brief Fix points for all parts. */
	std::vector<idx_t> _fixPoints;

	/** @brief Flags used to recognize whether the specified property was computed. */
	std::vector<bool> _flags;

	/** @brief Size of an element with maximal number of nodes. */
	size_t _maxElementSize;
};


class SurfaceMesh: public Mesh
{

public:
	SurfaceMesh(): Mesh() { };
	SurfaceMesh(const char *mesh, const char *coordinates, idx_t parts, idx_t fixPoints):
		Mesh(mesh, coordinates, parts, fixPoints) { };
	SurfaceMesh(const Mesh &mesh)
	{
		mesh.getSurface(*this);
	}

	void elasticity(DenseMatrix &K, size_t part) const;
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



#endif /* MESH_H_ */
