#ifndef TETRAHEDRON10_H_
#define TETRAHEDRON10_H_

#include "../element.h"
#include "../2D/triangle3.h"
#include "../2D/triangle6.h"
#include "tetrahedron4.h"

#define Tetrahedron10NodesCount 10
#define Tetrahedron10FacesCount 4
#define Tetrahedron10GPCount 15
#define Tetrahedron10VTKCode 24

namespace espreso {

class Tetrahedron10: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Tetrahedron10(const eslocal *indices, eslocal n, const eslocal *params);
	Tetrahedron10(std::ifstream &is);

	Element* copy() const
	{
		return new Tetrahedron10(*this);
	}

	eslocal vtkCode() const
	{
		return Tetrahedron10VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Tetrahedron10NodesCount;
	}

	size_t coarseSize() const
	{
		return Tetrahedron4NodesCount;
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

	const std::vector<Property>& DOFElement() const
	{
		return Tetrahedron10::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Tetrahedron10::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Tetrahedron10::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Tetrahedron10::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Tetrahedron10::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getCoarseFace(size_t face) const
	{
		return Tetrahedron4::getF(_indices, _params, face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Tetrahedron10NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;
};

}

#endif /* TETRAHEDRON10_H_ */
