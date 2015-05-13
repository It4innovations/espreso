#ifndef TETRAHEDRON_H_
#define TETRAHEDRON_H_

#include "../element.h"
#include "../2D/triangle.h"

#define TetrahedronNodesCount 4
#define TetrahedronFacesCount 4
#define TetrahedronGPCount 4
#define TetrahedronVTKCode 10

class Tetrahedron: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Tetrahedron(idx_t *indices);

	Element* copy() const
	{
		return new Tetrahedron(*this);
	}

	int vtkCode() const
	{
		return TetrahedronVTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return TetrahedronNodesCount;
	}

	size_t gpSize() const
	{
		return TetrahedronGPCount;
	}

	size_t faces() const
	{
		return TetrahedronFacesCount;
	}

	const std::vector<std::vector<double> >& dN() const
	{
		return Tetrahedron::_dN;
	}

	const std::vector<std::vector<double> >&  N() const
	{
		return Tetrahedron::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Tetrahedron::_weighFactor;
	}

	std::vector<idx_t> getNeighbours(size_t nodeIndex) const;
	std::vector<idx_t> getFace(size_t face) const;

protected:

	idx_t* indices()
	{
		return _indices;
	}

private:
	idx_t _indices[TetrahedronNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};

#endif /* TETRAHEDRON_H_ */
