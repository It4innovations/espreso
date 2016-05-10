
#ifndef PYRAMID5_H_
#define PYRAMID5_H_

#include "../element.h"
#include "../2D/triangle3.h"
#include "../2D/square4.h"

#define Pyramid5NodesCount 5
#define Pyramid5FacesCount 5
#define Pyramid5GPCount 8
#define Pyramid5VTKCode 14

namespace espreso {

class Pyramid5: public Element
{
	friend class Pyramid13;

public:
	static bool match(const eslocal *indices, eslocal n);

	Pyramid5(const eslocal *indices, eslocal n, const eslocal *params);
	Pyramid5(std::ifstream &is);


	Element* copy() const
	{
		return new Pyramid5(*this);
	}

	eslocal vtkCode() const
	{
		return Pyramid5VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Pyramid5NodesCount;
	}

	size_t coarseSize() const
	{
		return Pyramid5NodesCount;
	}

	size_t gpSize() const
	{
		return Pyramid5GPCount;
	}

	size_t faces() const
	{
		return Pyramid5FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Pyramid5::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Pyramid5::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Pyramid5::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getFullFace(size_t face) const
	{
		return getF(_indices, _params, face);
	}

	Element* getCoarseFace(size_t face) const
	{
		return getFullFace(face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;

protected:
	static Element* getF(const eslocal *indices, const eslocal *params, size_t face);

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Pyramid5NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* PYRAMID5_H_ */
