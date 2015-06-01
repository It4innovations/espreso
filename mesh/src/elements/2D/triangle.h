#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../element.h"
#include "../1D/line.h"

#define TriangleNodesCount 3
#define TriangleFacesCount 1
#define TriangleGPCount 3
#define TriangleVTKCode 5

namespace mesh {

class Triangle: public Element
{

public:
	static bool match(esint *indices, esint n);

	Triangle(esint *indices);

	Element* copy() const
	{
		return new Triangle(*this);
	}

	esint vtkCode() const
	{
		return TriangleVTKCode;
	}

	const esint* indices() const
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

	std::vector<esint> getNeighbours(size_t nodeIndex) const;
	std::vector<esint> getFace(size_t face) const;

protected:

	esint* indices()
	{
		return _indices;
	}

private:
	esint _indices[TriangleNodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}



#endif /* TRIANGLE_H_ */
