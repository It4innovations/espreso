#include "triangle.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Triangle::_dN;
std::vector<DenseMatrix> Triangle::_N;
std::vector<double> Triangle::_weighFactor;

bool Triangle::match(eslocal *indices, eslocal n)
{
	if (n != 3) {
		return false;
	}

	for (eslocal i = 0; i < TriangleNodesCount - 1; i++) {
		for (eslocal j = i + 1; j < TriangleNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Triangle::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	result.reserve(2);

	for (eslocal i = 0; i < TriangleNodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}

	return result;
}

std::vector<eslocal> Triangle::getFace(size_t face) const
{
	if (face < 2) {
		return std::vector<eslocal> (_indices + face, _indices + face + 2);
	} else {
		return std::vector<eslocal> ({_indices[2], _indices[0]});
	}
}

Triangle::Triangle(eslocal *indices, eslocal *params): Element(params)
{
	memcpy(_indices, indices, TriangleNodesCount * sizeof(eslocal));
}

