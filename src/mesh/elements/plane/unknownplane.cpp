
#include "unknownplane.h"
#include "../line/unknownline.h"

using namespace espreso;

void UnknownPlane::fillEdges()
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

	_edgeNodes.resize(domains.size());
	for (size_t i = 0; i < _indices.size(); i++) {
		for (size_t d = 0; d < _nodes[_indices[i]]->domains().size(); d++) {
			_edgeNodes[std::find(domains.begin(), domains.end(), _nodes[_indices[i]]->domains()[d]) - domains.begin()].push_back(_indices[i]);
		}
	}

	std::sort(_edgeNodes.begin(), _edgeNodes.end(), [] (const std::vector<eslocal> &i, const std::vector<eslocal> &j) { return i.size() < j.size(); });

	for (size_t f = 0; f < _edgeNodes.size() - 1; f++) {
		_edges.push_back(new UnknownLine(_nodes, _edgeNodes[f], _stiffnessMatrix));
	}
}





