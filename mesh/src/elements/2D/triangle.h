#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../element.h"
#include "../1D/line.h"

#define TriangleNodesCount 3
#define TriangleFacesCount 1
#define TriangleGPCount 3
#define TriangleVTKCode 5

class Triangle: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Triangle(idx_t *indices);

	Element* copy() const
	{
		return new Triangle(*this);
	}

	int vtkCode() const
	{
		return TriangleVTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return TriangleNodesCount;
	}

	size_t gpSize() const
	{
		return TriangleGPCount;
	}

	size_t faces() const
	{
		return TriangleFacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Triangle::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Triangle::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Triangle::_weighFactor;
	}

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[TriangleNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};



#endif /* TRIANGLE_H_ */
