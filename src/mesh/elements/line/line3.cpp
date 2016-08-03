
#include "line3.h"

using namespace espreso;

size_t Line3::_counter = 0;

std::vector<Property> Line3::_DOFElement;
std::vector<Property> Line3::_DOFFace;
std::vector<Property> Line3::_DOFEdge;
std::vector<Property> Line3::_DOFPoint;
std::vector<Property> Line3::_DOFMidPoint;

// TODO: Implement base functions
std::vector<DenseMatrix> Line3::_dN;
std::vector<DenseMatrix> Line3::_N;
std::vector<double> Line3::_weighFactor;

bool Line3::match(const eslocal *indices, eslocal n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line3::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

Line3::Line3(const eslocal *indices)
{
	memcpy(_indices, indices, Line3NodesCount * sizeof(eslocal));
}




