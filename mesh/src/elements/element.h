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
#include "../matrices/sparseDOKMatrix.h"

struct ElementCmp;
class Element;

typedef std::map<Element*, std::pair<int, int>, ElementCmp> BoundaryFaces;
typedef std::map<Element*, std::set<int> , ElementCmp> BoundaryLines;
typedef std::vector<std::set<int> > BoundaryNodes;

class Element
{
public:
	enum FLAGS {
		GLOBAL,
		LOCAL
	};

	inline static bool match(idx_t *indices, idx_t x, idx_t y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, Element &e);

	virtual ~Element() {};

	idx_t node(size_t index) const
	{
		return indices()[index];
	}

	idx_t localNode(size_t index) const
	{
		return localIndices()[index];
	}

	void fillNodes(idx_t *nodes, FLAGS indexing) const;
	void setLocalIndices(idx_t *mapping);
	void fillBoundaries(BoundaryNodes &nodes, int part) const;
	void coordinatesToVector(std::vector<double> &vector, const Coordinates &coordinates) const;

	void elasticity(std::vector<double> &Ke, std::vector<double> &Me, std::vector<double> &fe,
		const Coordinates &coordinates, std::vector<double> &inertia, double ex, double mi) const
	{
		_elaticity(Ke, Me, fe, coordinates, inertia, ex, mi, true);
	}

	void elasticity(std::vector<double> &Ke, std::vector<double> &fe,
		const Coordinates &coordinates, std::vector<double> &inertia, double ex, double mi) const
	{
		std::vector<double> Me;
		_elaticity(Ke, Me, fe, coordinates, inertia, ex, mi, false);
	}

	void addLocalValues(SparseDOKMatrix &K, SparseDOKMatrix &M, std::vector<double> &f,
		const std::vector<double> &Ke, const std::vector<double> &Me, const std::vector<double> &fe, int offset) const
	{
		_addLocalValues(K, M, f, Ke, Me, fe, offset, true);
	}

	void addLocalValues(SparseDOKMatrix &K, std::vector<double> &f,
		const std::vector<double> &Ke, const std::vector<double> &fe, int offset) const
	{
		SparseDOKMatrix M(0, 0);
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
	virtual void fillNeighbour(BoundaryNodes &nodes, int indexing) const = 0;
	virtual size_t faces() const = 0;
	virtual void fillFaces(BoundaryFaces &faces, int part) const = 0;
	virtual void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const = 0;
	virtual void fillLines(BoundaryLines &lines, int parts[]) const = 0;
	virtual const idx_t* indices() const = 0;
	virtual idx_t* localIndices() = 0;
	virtual const idx_t* localIndices() const = 0;

protected:
	bool isOnBorder(const BoundaryNodes &nodes, const idx_t *positions, idx_t n) const;

	void _elaticity(
		std::vector<double> &Ke,
		std::vector<double> &Me,
		std::vector<double> &fe,
		const Coordinates &coordinates,
		std::vector<double> &inertia,
		double ex,
		double mi,
		bool dynamic
	) const;

	void _addLocalValues(
		SparseDOKMatrix &K,
		SparseDOKMatrix &M,
		std::vector<double> &f,
		const std::vector<double> &Ke,
		const std::vector<double> &Me,
		const std::vector<double> &fe,
		int offset,
		bool dynamic
	) const;
};

struct ElementCmp {
	bool operator()(const Element *e1, const Element *e2) const
	{
		if (e1->size() == e2->size()) {
			return memcmp(e1->indices(), e2->indices(), sizeof(idx_t) * e1->size()) < 0;
		} else {
			return e1->size() < e2->size();
		}
	}
};


#endif /* ELEMENT_H_ */
