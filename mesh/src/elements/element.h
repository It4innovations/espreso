#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <fstream>

#include "1D/point.h"
#include "metis.h"
#include "mkl.h"

#include "../structures/coordinates.h"
#include "../matrices/matrices.h"

class Element;

typedef std::vector<std::set<int> > BoundaryNodes;

class Element
{
public:
	enum IndicesType {
		GLOBAL,
		LOCAL
	};

	inline static bool match(idx_t *indices, idx_t x, idx_t y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, const Element &e);
	friend void operator<<(double *nodeArray, const Element &e)
	{
		for (size_t i = 0; i < e.size(); i++) {
			nodeArray[i] = e.node(i);
		}
	}

	virtual ~Element() {};

	idx_t node(size_t index) const
	{
		return indices()[index];
	}

	void fillNodes(idx_t *nodes) const;
	void setLocalIndices(std::vector<idx_t> &mapping);
	void fillBoundaries(BoundaryNodes &nodes, int part) const;
	void coordinatesToVector(std::vector<double> &vector, const Coordinates &coordinates,
			IndicesType indicesType, size_t part) const;

	void elasticity(std::vector<double> &Ke, std::vector<double> &Me, std::vector<double> &fe,
		std::vector<double> &coordinates, std::vector<double> &inertia, double ex, double mi) const
	{
		_elaticity(Ke, Me, fe, coordinates, inertia, ex, mi, true);
	}

	void elasticity(std::vector<double> &Ke, std::vector<double> &fe,
			std::vector<double> &coordinates, std::vector<double> &inertia, double ex, double mi) const
	{
		std::vector<double> Me;
		_elaticity(Ke, Me, fe, coordinates, inertia, ex, mi, false);
	}

	void addLocalValues(SparseVVPMatrix &K, SparseVVPMatrix &M, std::vector<double> &f,
		const std::vector<double> &Ke, const std::vector<double> &Me, const std::vector<double> &fe, int offset) const
	{
		_addLocalValues(K, M, f, Ke, Me, fe, offset, true);
	}

	void addLocalValues(SparseVVPMatrix &K, std::vector<double> &f,
		const std::vector<double> &Ke, const std::vector<double> &fe, int offset) const
	{
		SparseVVPMatrix M(0, 0);
		std::vector<double> Me;
		_addLocalValues(K, M, f, Ke, Me, fe, offset, false);
	}

	virtual Element* copy() const = 0;

	virtual const std::vector<std::vector<double> >& dN() const = 0;
	virtual const std::vector<std::vector<double> >&  N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	// Virtual methods
	virtual int vtkCode() const = 0;
	virtual size_t size() const = 0;
	virtual size_t gpSize() const = 0;
	virtual size_t faces() const = 0;
	virtual std::vector<idx_t> getFace(size_t face) const = 0;
	virtual std::vector<idx_t> getNeighbours(size_t nodeIndex) const = 0;
	virtual const idx_t* indices() const = 0;

protected:
	virtual idx_t* indices() = 0;

	void _elaticity(
		std::vector<double> &Ke,
		std::vector<double> &Me,
		std::vector<double> &fe,
		std::vector<double> &coordinates,
		std::vector<double> &inertia,
		double ex,
		double mi,
		bool dynamic
	) const;

	void _addLocalValues(
		SparseVVPMatrix &K,
		SparseVVPMatrix &M,
		std::vector<double> &f,
		const std::vector<double> &Ke,
		const std::vector<double> &Me,
		const std::vector<double> &fe,
		int offset,
		bool dynamic
	) const;
};




#endif /* ELEMENT_H_ */
