#ifndef HEXAHEDRON8_H_
#define HEXAHEDRON8_H_

#include "../element.h"
#include "../2D/square.h"
#include "../1D/line.h"

#define Hexahedron8NodesCount 8
#define Hexahedron8FacesCount 6
#define Hexahedron8GPCount 8
#define Hexahedron8VTKCode 12

namespace mesh {

class Hexahedron8: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);

	Hexahedron8(eslocal *indices);
	Hexahedron8(std::ifstream &is);

	Element* copy() const
	{
		return new Hexahedron8(*this);
	}

	eslocal vtkCode() const
	{
		return Hexahedron8VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Hexahedron8NodesCount;
	}

	size_t gpSize() const
	{
		return Hexahedron8GPCount;
	}

	size_t faces() const
	{
		return Hexahedron8FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Hexahedron8::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Hexahedron8::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Hexahedron8::_weighFactor;
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
	inline void setFaceNodes(eslocal nodes[], eslocal face) const;

	eslocal _indices[Hexahedron8NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}

#endif /* HEXAHEDRON8_H_ */
