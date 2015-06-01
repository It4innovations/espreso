#include "square.h"

using namespace mesh;

// TODO: Implement base functions
std::vector<DenseMatrix> Square::_dN;
std::vector<DenseMatrix> Square::_N;
std::vector<double> Square::_weighFactor;

bool Square::match(esint *indices, esint n)
{
	if (n != 4) {
		return false;
	}

	for (esint i = 0; i < SquareNodesCount - 1; i++) {
		for (esint j = i + 1; j < SquareNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<esint> Square::getNeighbours(size_t nodeIndex) const
{
	std::vector<esint> result(2);

	result[0] = _indices[(nodeIndex + 1) % 4];
	result[1] = _indices[(nodeIndex + 3) % 4];

	return result;
}

std::vector<esint> Square::getFace(size_t face) const
{
	return std::vector<esint> (_indices, _indices + 4);
}

Square::Square(esint *indices)
{
	memcpy(_indices, indices, SquareNodesCount * sizeof(esint));
}



