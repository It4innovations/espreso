
#ifndef HEXAHEDRON20_H_
#define HEXAHEDRON20_H_

#include "../element.h"
#include "../2D/square.h"
#include "../1D/line.h"

#define Hexahedron20NodesCount 8
#define Hexahedron20FacesCount 6
#define Hexahedron20GPCount 8
#define Hexahedron20VTKCode 12

class Hexahedron20: public Element
{

public:
	static bool match(idx_t *indices, idx_t n);

	Hexahedron20(idx_t *indices);

	Element* copy() const
	{
		return new Hexahedron20(*this);
	}

	int vtkCode() const
	{
		return Hexahedron20VTKCode;
	}

	const idx_t* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Hexahedron20NodesCount;
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

  const std::vector< std::vector< double > >& rst()
  {
    return Hexahedron20::_rst;
  }
	const std::vector<DenseMatrix>&  N() const
	{
		return Hexahedron20::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Hexahedron20::_weighFactor;
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

	idx_t _indices[Hexahedron20NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

#endif /* HEXAHEDRON20_H_ */
