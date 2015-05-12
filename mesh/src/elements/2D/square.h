#ifndef SQUARE_H_
#define SQUARE_H_

#include "../element.h"
#include "../1D/line.h"

#define SquareNodesCount 4
#define SquareFacesCount 1
#define SquareGPCount 4
#define SquareVTKCode 9

class Square: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Square(idx_t *indices);

	Element* copy() const
	{
		return new Square(*this);
	}

	int vtkCode() const
	{
		return SquareVTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
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

	void fillNeighbour(BoundaryNodes &nodes) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[SquareNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};




#endif /* SQUARE_H_ */
