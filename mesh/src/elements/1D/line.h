#ifndef LINE_H_
#define LINE_H_


#include "../element.h"

#define LineNodesCount 2
#define LineFacesCount 0
#define LineGPCount 2
#define LineVTKCode 3

class Line: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Line(idx_t *indices);

	Element* copy() const
	{
		return new Line(*this);
	}

	int vtkCode() const
	{
		return LineVTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
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

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[LineNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};



#endif /* LINE_H_ */
