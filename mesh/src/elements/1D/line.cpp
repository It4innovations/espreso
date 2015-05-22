#include "line.h"

// TODO: Implement base functions
std::vector<DenseMatrix> Line::_dN;
std::vector<DenseMatrix> Line::_N;
std::vector<double> Line::_weighFactor;

bool Line::match(idx_t *indices, idx_t n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<idx_t> Line::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

std::vector<idx_t> Line::getFace(size_t face) const
{
	return std::vector<idx_t> ();
}

Line::Line(idx_t *indices)
{
	memcpy(_indices, indices, LineNodesCount * sizeof(idx_t));
}




