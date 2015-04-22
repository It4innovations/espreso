#ifndef SQUARE_H_
#define SQUARE_H_

#include "../element.h"
#include "../1D/line.h"

#define SquareNodesCount 4
#define SquareFacesCount 1
#define SquareGPCount 4

class Square: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Square(idx_t *indices);

	Element* copy() const
	{
		return new Square(*this);
	}

	const idx_t* indices() const
	{
		return _indices;
	}

	const idx_t* localIndices() const
	{
		return _localIndices;
	}

	idx_t* localIndices()
	{
		return _localIndices;
	}

	size_t size() const
	{
		return SquareNodesCount;
	}

	size_t gpSize() const
	{
		return SquareGPCount;
	}

	size_t faces() const
	{
		return SquareFacesCount;
	}

	const std::vector<std::vector<double> >& dN() const
	{
		return Square::_dN;
	}

	const std::vector<std::vector<double> >&  N() const
	{
		return Square::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Square::_weighFactor;
	}

	void fillNeighbour(BoundaryNodes &nodes, int indexing) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

private:
	idx_t _indices[SquareNodesCount];
	idx_t _localIndices[SquareNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};




#endif /* SQUARE_H_ */
