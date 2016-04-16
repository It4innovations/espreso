#ifndef TRIANGLE6_H_
#define TRIANGLE6_H_

#include "../1D/line2.h"
#include "../element.h"
#include "triangle3.h"

#define Triangle6NodesCount 6
#define Triangle6FacesCount 3
#define Triangle6GPCount 3
#define Triangle6VTKCode 22

namespace espreso {

class Triangle6: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);

	Triangle6(eslocal *indices, eslocal *params);

	Element* copy() const
	{
		return new Triangle6(*this);
	}

	eslocal vtkCode() const
	{
		return Triangle6VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Triangle6NodesCount;
	}

	size_t coarseSize() const
	{
		return Triangle3NodesCount;
	}

	size_t gpSize() const
	{
		return Triangle6GPCount;
	}

	size_t faces() const
	{
		return Triangle6FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Triangle6::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Triangle6::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Triangle6::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 2;
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
	eslocal _indices[Triangle6NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}



#endif /* TRIANGLE6_H_ */
