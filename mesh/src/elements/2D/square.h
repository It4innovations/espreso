#ifndef SQUARE_H_
#define SQUARE_H_

#include "../element.h"
#include "../1D/line.h"

#define SquareNodesCount 4
#define SquareFacesCount 1
#define SquareGPCount 4
#define SquareVTKCode 9

namespace mesh {

class Square: public Element
{

public:
	static bool match(esint *indices, esint n);

	Square(esint *indices);

	Element* copy() const
	{
		return new Square(*this);
	}

	esint vtkCode() const
	{
		return SquareVTKCode;
	}

	const esint* indices() const
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

	std::vector<esint> getNeighbours(size_t nodeIndex) const;
	std::vector<esint> getFace(size_t face) const;

protected:

	esint* indices()
	{
		return _indices;
	}

private:
	esint _indices[SquareNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* SQUARE_H_ */
