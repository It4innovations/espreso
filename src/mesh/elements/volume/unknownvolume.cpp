
#include "unknownvolume.h"
#include "../plane/unknownplane.h"

using namespace espreso;

void UnknownVolume::fillFaces()
{
	std::vector<eslocal> domains;
	for (size_t i = 0; i < _indices.size(); i++) {
		domains.insert(domains.end(), _nodes[_indices[i]]->domains().begin(), _nodes[_indices[i]]->domains().end());
	}
	std::sort(domains.begin(), domains.end());
	Esutils::removeDuplicity(domains);
	if (domains.size() == 1) {
		return;
	}

	_faceNodes.resize(domains.size());
	for (size_t i = 0; i < _indices.size(); i++) {
		for (size_t d = 0; d < _nodes[_indices[i]]->domains().size(); d++) {
			_faceNodes[std::find(domains.begin(), domains.end(), _nodes[_indices[i]]->domains()[d]) - domains.begin()].push_back(_indices[i]);
		}
	}

	std::sort(_faceNodes.begin(), _faceNodes.end(), [&] (const std::vector<eslocal> &i, const std::vector<eslocal> &j) {
		if (i.size() > nCommon() && j.size() > nCommon()) {
			return i.size() < j.size();
		} else {
			return i.size() > j.size();
		}
	});

	for (size_t f = 0; f < _faceNodes.size() && _faceNodes[f].size() < _indices.size(); f++) {
		_faces.push_back(new UnknownPlane(_nodes, _faceNodes[f], _stiffnessMatrix));
		_faces.back()->parentElements().push_back(this);
	}
}




