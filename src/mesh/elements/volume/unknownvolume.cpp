
#include "unknownvolume.h"
#include "../plane/unknownplane.h"

using namespace espreso;

void UnknownVolume::fillFaces()
{
	std::vector<Element*> elements;
	for (size_t i = 0; i < _indices.size(); i++) {
		elements.insert(elements.end(), _nodes[_indices[i]]->parentElements().begin(), _nodes[_indices[i]]->parentElements().end());
	}
	std::sort(elements.begin(), elements.end());
	Esutils::removeDuplicity(elements);

	for (size_t e = 0; e < elements.size(); e++) {
		if (elements[e]->domains()[0] != domains()[0]) {
			std::vector<eslocal> intersection(_indices.size());
			std::vector<eslocal> ind0(elements[e]->indices(), elements[e]->indices() + elements[e]->nodes());
			std::vector<eslocal> ind1 = _indices;
			std::sort(ind0.begin(), ind0.end());
			std::sort(ind1.begin(), ind1.end());
			auto it = std::set_intersection(ind0.begin(), ind0.end(), ind1.begin(), ind1.end(), intersection.begin());
			if (it - intersection.begin() > nCommon()) {
				_faceNodes.push_back(std::vector<eslocal>(intersection.begin(), it));
			}
		}
	}

	for (size_t i = 0; i < _faceNodes.size(); i++) {
		_faces.push_back(new UnknownPlane(_nodes, _faceNodes[i], _stiffnessMatrix));
		_faces.back()->parentElements().push_back(this);
	}
}




