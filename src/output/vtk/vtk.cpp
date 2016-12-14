
#include "vtk.h"

namespace espreso {
namespace store{


void VTK::computeCenters()
{
	_sCenters.resize(_mesh.parts());
	for (size_t p = 0; p < _mesh.coordinates().parts(); p++) {
		for (size_t i = 0; i < _mesh.coordinates().localSize(p); i++) {
			_sCenters[p] += _mesh.coordinates().get(i, p);
		}
		_sCenters[p] /= _mesh.coordinates().localSize(p);
	}

	for (size_t i = 0; i < _mesh.coordinates().clusterSize(); i++) {
		_cCenter += _mesh.coordinates()[i];
	}
	_cCenter /= _mesh.coordinates().clusterSize();
}

Point VTK::shrink(const Point &p, size_t part) const
{
	Point x = p;
	x = _sCenters[part] + (x - _sCenters[part]) * _output.domain_shrink_ratio;
	x = _cCenter + (x - _cCenter) * _output.cluster_shrink_ratio;
	return x;
}

Point VTK::shrink(const Point &p, const Point &sCenter, const Point &cCenter) const
{
	Point x = p;
	x = sCenter + (x - sCenter) * _output.domain_shrink_ratio;
	x = cCenter + (x - cCenter) * _output.cluster_shrink_ratio;
	return x;
}


}
}


