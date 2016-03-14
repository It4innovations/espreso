#ifndef SQUARE_H_
#define SQUARE_H_

#include "../element.h"
#include "../1D/line.h"

#define SquareNodesCount 4
#define SquareFacesCount 4
#define SquareGPCount 4
#define SquareVTKCode 9

namespace espreso {

class Square: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);

	Square(eslocal *indices, eslocal *params);

	Element* copy() const
	{
		return new Square(*this);
	}

	eslocal vtkCode() const
	{
		return SquareVTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return SquareNodesCount;
	}

	size_t coarseSize() const
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

	const std::vector<DenseMatrix>& dN() const
	{
		return Square::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Square::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Square::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 2;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[SquareNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* SQUARE_H_ */
