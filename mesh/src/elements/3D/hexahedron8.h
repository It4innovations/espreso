#ifndef HEXAHEDRON8_H_
#define HEXAHEDRON8_H_

#include "../element.h"
#include "../2D/square.h"
#include "../1D/line.h"

#define Hexahedron8NodesCount 8
#define Hexahedron8FacesCount 6
#define Hexahedron8GPCount 8
#define Hexahedron8VTKCode 12

class Hexahedron8: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Hexahedron8(idx_t *indices);

	Element* copy() const
	{
		return new Hexahedron8(*this);
	}

	int vtkCode() const
	{
		return Hexahedron8VTKCode;
	}

	const idx_t* indices() const
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

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	inline void setFaceNodes(idx_t nodes[], idx_t face) const;

	idx_t _indices[Hexahedron8NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

#endif /* HEXAHEDRON8_H_ */
