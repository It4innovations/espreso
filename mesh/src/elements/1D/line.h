#ifndef LINE_H_
#define LINE_H_


#include "../element.h"

#define LineNodesCount 2
#define LineFacesCount 0
#define LineGPCount 2

class Line: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Line(idx_t *indices);

	Element* copy() const
	{
		return new Line(*this);
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
		return LineNodesCount;
	}

	size_t gpSize() const
	{
		return LineGPCount;
	}

	size_t faces() const
	{
		return LineFacesCount;
	}

	const std::vector<std::vector<double> >& dN() const
	{
		return Line::_dN;
	}

	const std::vector<std::vector<double> >&  N() const
	{
		return Line::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Line::_weighFactor;
	}

	void fillNeighbour(BoundaryNodes &nodes, int indexing) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

private:
	idx_t _indices[LineNodesCount];
	idx_t _localIndices[LineNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};



#endif /* LINE_H_ */
