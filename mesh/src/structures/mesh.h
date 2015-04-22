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

namespace flags
{

enum FLAGS {
	PARTITIONS,
	FIX_POINTS,
	BOUNDARIES,
	FLAGS_SIZE
};
}

typedef std::pair<int, double> IV_elem;

class Mesh
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);

	Mesh(const char *fileName, Coordinates &coordinates, idx_t parts, idx_t fixPoints);
	Mesh(size_t size, Coordinates &coordinates);

	void push_element(Element* e)
	{
		_elements.push_back(e);
	}

	~Mesh();

	void partitiate(idx_t parts, idx_t fixPoints);
	void computeFixPoints(idx_t fixPoints);

	const std::vector<Element*>& getElements() const
	{
		return _elements;
	};

	const std::vector<idx_t>& getPartition() const
	{
		return _partPtrs;
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
		SparseDOKMatrix _K(K.type(), K.rows(), K.columns());
		SparseDOKMatrix _M(M.type(), M.rows(), M.columns());
		_assemble_matrix(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}

	void assemble_matrix(SparseIJVMatrix &K, SparseIJVMatrix &M, std::vector<double> &f, idx_t part)
	{
		SparseDOKMatrix _K(K.type(), K.rows(), K.columns());
		SparseDOKMatrix _M(M.type(), M.rows(), M.columns());
		_assemble_matrix(_K, _M, f, part, true);
		K = _K;
		M = _M;
	}

	void assemble_matrix(SparseDOKMatrix &K, SparseDOKMatrix &M, std::vector<double> &f, idx_t part)
	{
		_assemble_matrix(K, M, f, part, true);
	}

	void assemble_matrix(SparseCSRMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseDOKMatrix _K(K.type(), K.rows(), K.columns());
		SparseDOKMatrix _M(0, 0);
		_assemble_matrix(_K, _M, f, part, false);
		K = _K;
	}

	void assemble_matrix(SparseIJVMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseDOKMatrix _K(K.type(), K.rows(), K.columns());
		SparseDOKMatrix _M(0, 0);
		_assemble_matrix(_K, _M, f, part, false);
		K = _K;
	}

	void assemble_matrix(SparseDOKMatrix &K, std::vector<double> &f, idx_t part)
	{
		SparseDOKMatrix M(0, 0);
		_assemble_matrix(K, M, f, part, false);
	}

private:
	Element* createElement(idx_t *indices, idx_t n);

	void _assemble_matrix(SparseDOKMatrix &K, SparseDOKMatrix &M, std::vector<double> &f, idx_t part, bool dynamic);

	idx_t* getPartition(idx_t first, idx_t last, idx_t parts) const;
	idx_t getCentralNode(idx_t first, idx_t last, idx_t *ePartition, idx_t part, idx_t subpart) const;

	void partitiate(idx_t *ePartition);
	void computeLocalIndices();

	void checkMETISResult(int result) const;
	void checkMKLResult(MKL_INT result) const;

	/** @brief Reference to coordinates. */
	Coordinates &_coordinates;

	/** @brief Array that stores all elements of the mesh. */
	std::vector<Element*> _elements;

	/** @brief The biggest node's id in the mesh. */
	idx_t _lastNode;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<idx_t> _partPtrs;

	/** @brief Number of nodes in each part. */
	std::vector<idx_t> _partsNodesCount;

	/** @brief Keeps mapping of nodes to mesh parts. */
	BoundaryNodes _boundaries;

	/** @brief Fix points for all parts. */
	std::vector<idx_t> _fixPoints;

	/** @brief Flags used to recognize whether the specified property was computed. */
	std::vector<bool> _flags;

	/** @brief Size of an element with maximal number of nodes. */
	size_t _maxElementSize;
};



#endif /* MESH_H_ */
