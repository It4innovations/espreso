#include "triangle.h"

using namespace mesh;

// TODO: Implement base functions
std::vector<DenseMatrix> Triangle::_dN;
std::vector<DenseMatrix> Triangle::_N;
std::vector<double> Triangle::_weighFactor;

bool Triangle::match(esint *indices, esint n)
{
	if (n != 3) {
		return false;
	}

	for (esint i = 0; i < TriangleNodesCount - 1; i++) {
		for (esint j = i + 1; j < TriangleNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<esint> Triangle::getNeighbours(size_t nodeIndex) const
{
	std::vector<esint> result;
	result.reserve(2);

	for (esint i = 0; i < TriangleNodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}

	return result;
}

std::vector<esint> Triangle::getFace(size_t face) const
{
	return std::vector<esint> (_indices, _indices + 3);
}

Triangle::Triangle(esint *indices)
{
	memcpy(_indices, indices, TriangleNodesCount * sizeof(esint));
}

