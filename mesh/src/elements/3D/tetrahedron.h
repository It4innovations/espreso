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

	const idx_t* localIndices() const
	{
		return _localIndices;
	}

	idx_t* localIndices()
	{
		return _localIndices;
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

	void fillNeighbour(BoundaryNodes &nodes, int indexing) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

private:
	inline void setFaceNodes(idx_t nodes[], idx_t face) const;

	idx_t _indices[TetrahedronNodesCount];
	idx_t _localIndices[TetrahedronNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};

#endif /* TETRAHEDRON_H_ */
