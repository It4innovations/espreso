#ifndef TETRAHEDRON10_H_
#define TETRAHEDRON10_H_

#include "../element.h"
#include "../2D/triangle.h"

#define Tetrahedron10NodesCount 10
#define Tetrahedron10FacesCount 4
#define Tetrahedron10GPCount 15
#define Tetrahedron10VTKCode 24

namespace mesh {

class Tetrahedron10: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Tetrahedron10(idx_t *indices);

	Element* copy() const
	{
		return new Tetrahedron10(*this);
	}

	int vtkCode() const
	{
		return Tetrahedron10VTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Tetrahedron10NodesCount;
	}

	size_t gpSize() const
	{
		return Tetrahedron10GPCount;
	}

	size_t faces() const
	{
		return Tetrahedron10FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Tetrahedron10::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Tetrahedron10::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Tetrahedron10::_weighFactor;
	}

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[Tetrahedron10NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}

#endif /* TETRAHEDRON10_H_ */
