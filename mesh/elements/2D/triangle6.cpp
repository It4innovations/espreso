#include "triangle6.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Triangle6::_dN;
std::vector<DenseMatrix> Triangle6::_N;
std::vector<double> Triangle6::_weighFactor;

bool Triangle6::match(eslocal *indices, eslocal n)
{
	if (n != 6) {
		return false;
	}

	for (eslocal i = 0; i < Triangle6NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Triangle6NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Triangle6::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	if (nodeIndex < 3) {
		result[0] = _indices[nodeIndex + 3];
		result[1] = _indices[(nodeIndex + 2) % 3 + 3];
	} else {
		result[0] = _indices[(nodeIndex + 1) % 3];
		result[1] = _indices[nodeIndex - 3];
	}

	return result;
}

std::vector<eslocal> Triangle6::getFace(size_t face) const
{
	if (face < 2) {
		return std::vector<eslocal> (_indices + face, _indices + face + 2);
	} else {
		return std::vector<eslocal> ({_indices[2], _indices[0]});
	}
}

Triangle6::Triangle6(eslocal *indices, eslocal *params): Element(params)
{
	memcpy(_indices, indices, Triangle6NodesCount * sizeof(eslocal));
}

