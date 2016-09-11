
#include "vtk.h"

using namespace espreso::output;

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

espreso::Point VTK::shrink(const espreso::Point &p, size_t part)
{
	espreso::Point x = p;
	x = _sCenters[part] + (x - _sCenters[part]) * _shrinkSubdomain;
	x = _cCenter + (x - _cCenter) * _shringCluster;
	return x;
}




