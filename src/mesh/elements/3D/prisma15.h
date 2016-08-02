
#ifndef PRISMA15_H_
#define PRISMA15_H_

#include "../element.h"
#include "prisma6.h"
#include "../2D/triangle6.h"
#include "../2D/square8.h"

#define Prisma15NodesCount 15
#define Prisma15FacesCount 5
#define Prisma15GPCount 9
#define Prisma15VTKCode 26

namespace espreso {

class Prisma15: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Prisma15(const eslocal *indices, eslocal n, const eslocal *params);
	Prisma15(std::ifstream &is);

	Element* copy() const
	{
		return new Prisma15(*this);
	}

	eslocal vtkCode() const
	{
		return Prisma15VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Prisma15NodesCount;
	}

	size_t coarseSize() const
	{
		return Prisma6NodesCount;
	}

	size_t gpSize() const
	{
		return Prisma15GPCount;
	}

	size_t faces() const
	{
		return Prisma15FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Prisma15::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Prisma15::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Prisma15::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Prisma15::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Prisma15::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Prisma15::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Prisma15::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Prisma15::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getCoarseFace(size_t face) const
	{
		return Prisma6::getF(_indices, _params, face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Prisma15NodesCount];

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
