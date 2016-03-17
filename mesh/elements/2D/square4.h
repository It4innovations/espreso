#ifndef SQUARE4_H_
#define SQUARE4_H_

#include "../element.h"
#include "../1D/line.h"

#define Square4NodesCount 4
#define Square4FacesCount 4
#define Square4GPCount 4
#define Square4VTKCode 9

namespace espreso {

class Square4: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Square4(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Square4(*this);
	}

	eslocal vtkCode() const
	{
		return Square4VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Square4NodesCount;
	}

	size_t coarseSize() const
	{
		return Square4NodesCount;
	}

	size_t gpSize() const
	{
		return Square4GPCount;
	}

	size_t faces() const
	{
		return Square4FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Square4::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Square4::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Square4::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 2;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;
	Element* getCoarseFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Square4NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* SQUARE4_H_ */
