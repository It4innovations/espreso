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

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[SquareNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};




#endif /* SQUARE_H_ */
