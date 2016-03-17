#ifndef LINE_H_
#define LINE_H_


#include "../element.h"

#define LineNodesCount 2
#define LineFacesCount 0
#define LineGPCount 2
#define LineVTKCode 3

namespace espreso {

class Line: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);

	Line(eslocal *indices, eslocal *params);

	Element* copy() const
	{
		return new Line(*this);
	}

	eslocal vtkCode() const
	{
		return LineVTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return LineNodesCount;
	}

	size_t coarseSize() const
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

	const std::vector<DenseMatrix>& dN() const
	{
		return Line::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Line::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Line::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 1;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;
	Element* getCoarseFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[LineNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* LINE_H_ */
