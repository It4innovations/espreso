#include "boundaries.h"

using namespace mesh;

Boundaries::Boundaries(const Mesh &m): _boundaries(m.coordinates().size()), _mesh(m)
{
	const std::vector<eslocal> &parts = m.getPartition();
	const std::vector<Element*> &elements = m.getElements();
	const Coordinates &c = m.coordinates();

	for (size_t p = 0; p + 1 < parts.size(); p++) {
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			for (size_t n = 0; n < elements[e]->size(); n++) {
				_boundaries[c.clusterIndex(elements[e]->node(n), p)].insert(p);
			}
		}
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Boundaries &b)
{
	std::set<eslocal>::const_iterator it;

	for (size_t i = 0; i < b._boundaries.size(); i++) {
		os << b._mesh.coordinates()[i] << ": ";
		for (it = b._boundaries[i].begin(); it != b._boundaries[i].end(); ++it) {
			os << *it << " ";
		}
		os << "\n";
	}
	return os;
}



