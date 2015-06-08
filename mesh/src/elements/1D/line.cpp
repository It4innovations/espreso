#include "line.h"

using namespace mesh;

// TODO: Implement base functions
std::vector<DenseMatrix> Line::_dN;
std::vector<DenseMatrix> Line::_N;
std::vector<double> Line::_weighFactor;

bool Line::match(eslocal *indices, eslocal n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

std::vector<eslocal> Line::getFace(size_t face) const
{
	return std::vector<eslocal> ();
}

Line::Line(eslocal *indices)
{
	memcpy(_indices, indices, LineNodesCount * sizeof(eslocal));
}




