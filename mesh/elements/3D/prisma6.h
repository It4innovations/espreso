
#ifndef PRISMA6_H_
#define PRISMA6_H_

#include "../element.h"
#include "../1D/line.h"

#define Prisma6NodesCount 6
#define Prisma6FacesCount 5
#define Prisma6GPCount 9
#define Prisma6VTKCode 13

namespace espreso {

class Prisma6: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

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

	eslocal nCommon() const
	{
		return 3;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Prisma6NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* PRISMA6_H_ */
