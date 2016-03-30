#ifndef TRIANGLE3_H_
#define TRIANGLE3_H_

#include "../element.h"
#include "../1D/line.h"

#define Triangle3NodesCount 3
#define Triangle3FacesCount 3
#define Triangle3GPCount 3
#define Triangle3VTKCode 5

namespace espreso {

class Triangle3: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);

	Triangle3(eslocal *indices, eslocal *params);

	Element* copy() const
	{
		return new Triangle3(*this);
	}

	eslocal vtkCode() const
	{
		return Triangle3VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Triangle3NodesCount;
	}

	size_t coarseSize() const
	{
		return Triangle3NodesCount;
	}

	size_t gpSize() const
	{
		return Triangle3GPCount;
	}

	size_t faces() const
	{
		return Triangle3FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Triangle3::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Triangle3::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Triangle3::_weighFactor;
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
	eslocal _indices[Triangle3NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}



#endif /* TRIANGLE3_H_ */
