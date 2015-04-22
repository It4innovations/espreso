#ifndef HEXAHEDRON_H_
#define HEXAHEDRON_H_

#include "../element.h"
#include "../2D/square.h"
#include "../1D/line.h"

#define HexahedronNodesCount 8
#define HexahedronFacesCount 6
#define HexahedronGPCount 8

class Hexahedron: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Hexahedron(idx_t *indices);

	Element* copy() const
	{
		return new Hexahedron(*this);
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
		return HexahedronNodesCount;
	}

	size_t gpSize() const
	{
		return HexahedronGPCount;
	}

	size_t faces() const
	{
		return HexahedronFacesCount;
	}

	const std::vector<std::vector<double> >& dN() const
	{
		return Hexahedron::_dN;
	}

	const std::vector<std::vector<double> >&  N() const
	{
		return Hexahedron::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Hexahedron::_weighFactor;
	}

	void fillNeighbour(BoundaryNodes &nodes, int indexing) const;
	void fillFaces(BoundaryFaces &faces, int part) const;
	void fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const;
	void fillLines(BoundaryLines &lines, int parts[]) const;

private:
	inline void setFaceNodes(idx_t nodes[], idx_t face) const;

	idx_t _indices[HexahedronNodesCount];
	idx_t _localIndices[HexahedronNodesCount];

	static std::vector<std::vector<double> > _dN;
	static std::vector<std::vector<double> > _N;
	static std::vector<double> _weighFactor;
};

#endif /* HEXAHEDRON_H_ */
