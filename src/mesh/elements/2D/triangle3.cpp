#include "triangle3.h"

using namespace espreso;

// TODO: Implement base functions
std::vector<DenseMatrix> Triangle3::_dN;
std::vector<DenseMatrix> Triangle3::_N;
std::vector<double> Triangle3::_weighFactor;

bool Triangle3::match(eslocal *indices, eslocal n)
{
	if (n != 3) {
		return false;
	}

	for (eslocal i = 0; i < Triangle3NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Triangle3NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Triangle3::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	result.reserve(2);

	for (eslocal i = 0; i < Triangle3NodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}

	return result;
}

std::vector<eslocal> Triangle3::getFace(size_t face) const
{
	if (face < 2) {
		return std::vector<eslocal> (_indices + face, _indices + face + 2);
	} else {
		return std::vector<eslocal> ({_indices[2], _indices[0]});
	}
}

Element* Triangle3::getFullFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Element* Triangle3::getCoarseFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Triangle3::Triangle3(const eslocal *indices, const eslocal *params): Element(params)
{
	memcpy(_indices, indices, Triangle3NodesCount * sizeof(eslocal));
}

