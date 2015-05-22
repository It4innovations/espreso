#include "square.h"

// TODO: Implement base functions
std::vector<DenseMatrix> Square::_dN;
std::vector<DenseMatrix> Square::_N;
std::vector<double> Square::_weighFactor;

bool Square::match(idx_t *indices, idx_t n)
{
	if (n != 4) {
		return false;
	}

	for (idx_t i = 0; i < SquareNodesCount - 1; i++) {
		for (idx_t j = i + 1; j < SquareNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Square::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result(2);

	result[0] = _indices[(nodeIndex + 1) % 4];
	result[1] = _indices[(nodeIndex + 3) % 4];

	return result;
}

std::vector<idx_t> Square::getFace(size_t face) const
{
	return std::vector<idx_t> (_indices, _indices + 4);
}

Square::Square(idx_t *indices)
{
	memcpy(_indices, indices, SquareNodesCount * sizeof(idx_t));
}



