#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../element.h"
#include "../1D/line.h"

#define TriangleNodesCount 3
#define TriangleFacesCount 1
#define TriangleGPCount 3
#define TriangleVTKCode 5

class Triangle: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Triangle(idx_t *indices);

	Element* copy() const
	{
		return new Triangle(*this);
	}

	int vtkCode() const
	{
		return TriangleVTKCode;
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
		return TriangleNodesCount;
	}

	size_t gpSize() const
	{
		return TriangleGPCount;
	}

	size_t faces() const
	{
		return TriangleFacesCount;
	}

	const std::vector<std::vector<double> >& dN() const
	{
		return Triangle::_dN;
	}

	const std::vector<std::vector<double> >&  N() const
	{
		return Triangle::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Triangle::_weighFactor;
	}

	void fillNeighbour(BoundaryNodes &nodes, int indexing) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

private:
	idx_t _indices[TriangleNodesCount];
	idx_t _localIndices[TriangleNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};



#endif /* TRIANGLE_H_ */
