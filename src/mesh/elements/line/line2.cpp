
#include "line2.h"

using namespace espreso;

size_t Line2::_counter = 0;

std::vector<Property> Line2::_DOFElement;
std::vector<Property> Line2::_DOFFace;
std::vector<Property> Line2::_DOFEdge;
std::vector<Property> Line2::_DOFPoint;
std::vector<Property> Line2::_DOFMidPoint;

// TODO: Implement base functions
std::vector<DenseMatrix> Line2::_dN;
std::vector<DenseMatrix> Line2::_N;
std::vector<double> Line2::_weighFactor;

bool Line2::match(const eslocal *indices, eslocal n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line2::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

Line2::Line2(const eslocal *indices)
{
	memcpy(_indices, indices, Line2NodesCount * sizeof(eslocal));
}




