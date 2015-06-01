
#ifndef TETRAHEDRON4_H_
#define TETRAHEDRON4_H_

#include "../element.h"
#include "../2D/triangle.h"

#define Tetrahedron4NodesCount 4
#define Tetrahedron4FacesCount 4
#define Tetrahedron4GPCount 4
#define Tetrahedron4VTKCode 10

namespace mesh {

class Tetrahedron4: public Element
{

public:
	static bool match(esint *indices, esint n);

	Tetrahedron4(esint *indices);

	Element* copy() const
	{
		return new Tetrahedron4(*this);
	}

	esint vtkCode() const
	{
		return Tetrahedron4VTKCode;
	}

	const esint* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Tetrahedron4NodesCount;
	}

	size_t gpSize() const
	{
		return Tetrahedron4GPCount;
	}

	size_t faces() const
	{
		return Tetrahedron4FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Tetrahedron4::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Tetrahedron4::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Tetrahedron4::_weighFactor;
	}

	std::vector<esint> getNeighbours(size_t nodeIndex) const;
	std::vector<esint> getFace(size_t face) const;

protected:

	esint* indices()
	{
		return _indices;
	}

private:
	esint _indices[Tetrahedron4NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};


}

#endif /* TETRAHEDRON4_H_ */
