
#ifndef PRISMA6_H_
#define PRISMA6_H_

#include "../element.h"
#include "../plane/triangle3.h"
#include "../plane/square4.h"

#define Prisma6NodesCount 6
#define Prisma6FacesCount 5
#define Prisma6GPCount 9
#define Prisma6VTKCode 13

namespace espreso {

class Prisma6: public Element
{
	friend class Prisma15;

public:
	static bool match(const eslocal *indices, eslocal n);
	static size_t counter()
	{
		return _counter;
	}
	static void setDOFs(
			const std::vector<Property> element,
			const std::vector<Property> face,
			const std::vector<Property> edge,
			const std::vector<Property> point,
			const std::vector<Property> midPoint)
	{
		_DOFElement = element;
		_DOFFace = face;
		_DOFEdge = edge;
		_DOFPoint = point;
		_DOFMidPoint = midPoint;
	}

	Prisma6(const eslocal *indices, eslocal n, const eslocal *params);
	Prisma6(std::ifstream &is);

	Element* copy() const
	{
		return new Prisma6(*this);
	}

	eslocal vtkCode() const
	{
		return Prisma6VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Prisma6NodesCount;
	}

	size_t coarseSize() const
	{
		return Prisma6NodesCount;
	}

	size_t gpSize() const
	{
		return Prisma6GPCount;
	}

	size_t faces() const
	{
		return Prisma6FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Prisma6::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Prisma6::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Prisma6::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Prisma6::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Prisma6::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Prisma6::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Prisma6::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Prisma6::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getFullFace(size_t face) const
	{
		return getF(_indices, _params, face);
	}

	Element* getCoarseFace(size_t face) const
	{
		return getFullFace(face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;

protected:
	static Element* getF(const eslocal *indices, const eslocal *params, size_t face);

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Prisma6NodesCount];

	static size_t _counter;

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


#endif /* PRISMA6_H_ */
