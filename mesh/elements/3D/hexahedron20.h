
#ifndef HEXAHEDRON20_H_
#define HEXAHEDRON20_H_

#include "../element.h"
#include "../2D/square4.h"
#include "../1D/line.h"
#include "hexahedron8.h"

#define Hexahedron20NodesCount 20
#define Hexahedron20FacesCount 6
#define Hexahedron20GPCount 8
#define Hexahedron20VTKCode 25

namespace espreso {

class Hexahedron20: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Hexahedron20(const eslocal *indices, eslocal n, const eslocal *params);
	Hexahedron20(std::ifstream &is);

	Element* copy() const
	{
		return new Hexahedron20(*this);
	}

	eslocal vtkCode() const
	{
		return Hexahedron20VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Hexahedron20NodesCount;
	}

	size_t coarseSize() const
	{
		return Hexahedron8NodesCount;
	}

	size_t gpSize() const
	{
		return Hexahedron20GPCount;
	}

	size_t faces() const
	{
		return Hexahedron20FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Hexahedron20::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Hexahedron20::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Hexahedron20::_weighFactor;
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
	eslocal _indices[Hexahedron20NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};


}

#endif /* HEXAHEDRON20_H_ */
