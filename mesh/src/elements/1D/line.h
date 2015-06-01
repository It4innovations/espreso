#ifndef LINE_H_
#define LINE_H_


#include "../element.h"

#define LineNodesCount 2
#define LineFacesCount 0
#define LineGPCount 2
#define LineVTKCode 3

namespace mesh {

class Line: public Element
{

public:
	static bool match(esint *indices, esint n);

	Line(esint *indices);

	Element* copy() const
	{
		return new Line(*this);
	}

	esint vtkCode() const
	{
		return LineVTKCode;
	}

	const esint* indices() const
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

	std::vector<esint> getNeighbours(size_t nodeIndex) const;
	std::vector<esint> getFace(size_t face) const;

protected:

	esint* indices()
	{
		return _indices;
	}

private:
	esint _indices[LineNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* LINE_H_ */
