#include "square4.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Square4::_dN;
std::vector<DenseMatrix> Square4::_N;
std::vector<double> Square4::_weighFactor;

bool Square4::match(const eslocal *indices, eslocal n)
{
	if (n != 4) {
		return false;
	}

	for (eslocal i = 0; i < Square4NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Square4NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Square4::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	result[0] = _indices[(nodeIndex + 1) % 4];
	result[1] = _indices[(nodeIndex + 3) % 4];

	return result;
}

std::vector<eslocal> Square4::getFace(size_t face) const
{
	if (face < 3) {
		return std::vector<eslocal>(_indices + face, _indices + face + 2);
	} else {
		return std::vector<eslocal>({_indices[3], _indices[0]});
	}
}

Element* getF(const eslocal *indices, const eslocal *params, size_t face)
{
	std::vector<eslocal> result(2);
	if (face < 3) {
		result[0] = indices[face];
		result[1] = indices[face + 1];
	} else {
		result[0] = indices[3];
		result[1] = indices[0];
	}

	return new Line2(indices, params);
}

Element* Square4::getFullFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Element* Square4::getCoarseFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Square4::Square4(const eslocal *indices, const eslocal *params): Element(params)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
}



