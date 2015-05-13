#include "triangle.h"

// TODO: Implement base functions
std::vector<std::vector<double> > Triangle::_dN;
std::vector<std::vector<double> > Triangle::_N;
std::vector<double> Triangle::_weighFactor;

bool Triangle::match(idx_t *indices, idx_t n)
{
	if (n != 3) {
		return false;
	}

	for (idx_t i = 0; i < TriangleNodesCount - 1; i++) {
		for (idx_t j = i + 1; j < TriangleNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Triangle::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result;
	result.reserve(2);

	for (idx_t i = 0; i < TriangleNodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}

	return result;
}

std::vector<idx_t> Triangle::getFace(size_t face) const
{
	return std::vector<idx_t> (_indices, _indices + 3);
}

Triangle::Triangle(idx_t *indices)
{
	memcpy(_indices, indices, TriangleNodesCount * sizeof(idx_t));
}

