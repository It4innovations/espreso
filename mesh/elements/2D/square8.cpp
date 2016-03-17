#include "square8.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Square8::_dN;
std::vector<DenseMatrix> Square8::_N;
std::vector<double> Square8::_weighFactor;

bool Square8::match(eslocal *indices, eslocal n)
{
	if (n != 8) {
		return false;
	}

	for (eslocal i = 0; i < Square8NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Square8NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Square8::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	if (nodeIndex < 4) {
		result[0] = _indices[nodeIndex + 4];
		result[1] = _indices[(nodeIndex + 3) % 4 + 4];
	} else {
		result[0] = _indices[(nodeIndex + 5) % 4];
		result[1] = _indices[nodeIndex - 4];
	}

	return result;
}

std::vector<eslocal> Square8::getFace(size_t face) const
{
	std::vector<eslocal> result(3);
	if (face < 3) {
		result[0] = _indices[face];
		result[1] = _indices[face + 1];
		result[2] = _indices[face + 4];
	} else {
		result[0] = _indices[3];
		result[1] = _indices[0];
		result[2] = _indices[7];
	}
	return result;
}

Element* Square8::getFullFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Element* Square8::getCoarseFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Square8::Square8(eslocal *indices, eslocal *params): Element(params)
{
	memcpy(_indices, indices, Square8NodesCount * sizeof(eslocal));
}



