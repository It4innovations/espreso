#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>

#include "mkl.h"
#include "cilk/cilk.h"

#include "../elements/elements.h"
#include "coordinates.h"

#include "../matrices/sparseDOKMatrix.h"
#include "../settings.h"

#include <vector>

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

typedef std::pair<int, double> IV_elem;

class Mesh
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);

	Mesh(Coordinates &coordinates);
	Mesh(const char *fileName, Coordinates &coordinates, idx_t parts, idx_t fixPoints);

	Mesh(const Mesh &other);
	Mesh& operator=(const Mesh &other);

	void saveVTK(std::vector<std::vector <int> > &l2g_vec);
	void saveVTK(std::vector<std::vector<double> > &displacement, std::vector<std::vector <int> > &l2g_vec );

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

	void assemble_matrix(SparseCSRMatrix &K, SparseCSRMatrix &M, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_assemble_matrix(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}

	void assemble_matrix(SparseIJVMatrix &K, SparseIJVMatrix &M, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_assemble_matrix(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}

	void assemble_matrix(SparseCSRMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_assemble_matrix(_K, _M, f, part, false);
		K = _K;
	}

	void assemble_matrix(SparseIJVMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseVVPMatrix _K;
		SparseVVPMatrix _M;
		_assemble_matrix(_K, _M, f, part, false);
		K = _K;
	}

private:
	static void assign(Mesh &m1, Mesh &m2);

	void saveBasis(std::ofstream &vtk, std::vector<std::vector<int> > &l2g_vec);

	Element* createElement(idx_t *indices, idx_t n);

	void _assemble_matrix(SparseVVPMatrix &K, SparseVVPMatrix &M, std::vector<double> &f, idx_t part, bool dynamic);

	idx_t* getPartition(idx_t first, idx_t last, idx_t parts) const;
	idx_t getCentralNode(idx_t first, idx_t last, idx_t *ePartition, idx_t part, idx_t subpart) const;

	void partitiate(idx_t *ePartition);
	void computeLocalIndices(size_t part);

	void checkMETISResult(int result) const;
	void checkMKLResult(MKL_INT result) const;

	/** @brief Reference to coordinates. */
	Coordinates &_coordinates;

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



#endif /* MESH_H_ */
