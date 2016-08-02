#ifndef TRIANGLE6_H_
#define TRIANGLE6_H_

#include "../1D/line2.h"
#include "../element.h"
#include "triangle3.h"

#define Triangle6NodesCount 6
#define Triangle6FacesCount 3
#define Triangle6GPCount 6
#define Triangle6VTKCode 22

namespace espreso {

class Triangle6: public Element
{

public:
	static bool match(const eslocal *indices, const eslocal n);

	Triangle6(const eslocal *indices, const eslocal *params);

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

	const std::vector<Property>& DOFElement() const
	{
		return Triangle6::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Triangle6::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Triangle6::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Triangle6::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Triangle6::_DOFMidPoint;
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

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;
};

}



#endif /* TRIANGLE6_H_ */
