
#ifndef PYRAMID13_H_
#define PYRAMID13_H_

#include "../element.h"
#include "../1D/line.h"
#include "pyramid5.h"

#define Pyramid13NodesCount 13
#define Pyramid13FacesCount 5
#define Pyramid13GPCount 8
#define Pyramid13VTKCode 27

namespace espreso {

class Pyramid13: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Pyramid13(const eslocal *indices, eslocal n, const eslocal *params);
	Pyramid13(std::ifstream &is);

	Element* copy() const
	{
		return new Pyramid13(*this);
	}

	eslocal vtkCode() const
	{
		return Pyramid13VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Pyramid13NodesCount;
	}

	size_t coarseSize() const
	{
		return Pyramid5NodesCount;
	}

	size_t gpSize() const
	{
		return Pyramid13GPCount;
	}

	size_t faces() const
	{
		return Pyramid13FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Pyramid13::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Pyramid13::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Pyramid13::_weighFactor;
	}

	eslocal nCommon() const
	{
		return 3;
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
	eslocal _indices[Pyramid13NodesCount];

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;
};

}


#endif /* PYRAMID5_H_ */
